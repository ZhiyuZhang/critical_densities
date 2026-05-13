"""Calculate optically thin, no-background critical densities after Shirley (2015)."""  # Module purpose.

from dataclasses import dataclass  # Small immutable-like containers keep parsed LAMDA data readable.
from pathlib import Path  # Path is used so data files are found relative to this script.

import numpy as np  # NumPy provides arrays, interpolation, and vectorized arithmetic.


CM_TO_K = 1.438776877  # Convert energies from cm^-1 to K using h*c/k_B.
DATA_DIR = Path(__file__).resolve().parent  # Directory containing this script and the .dat files.
MOLECULE_FILES = {  # Map accepted molecule names to local LAMDA data files.
    "12CO": "12co.dat",  # Main CO isotopologue.
    "13CO": "13co.dat",  # 13CO isotopologue.
    "C18O": "c18o.dat",  # C18O isotopologue.
    "HCN": "hcn.dat",  # HCN data.
    "HCO+": "hco+.dat",  # HCO+ data.
    "HNC": "hnc.dat",  # HNC data.
    "CS": "cs.dat",  # CS data.
    "CN": "cn.dat",  # CN fine-structure data.
}  # End molecule map.


@dataclass(frozen=True)  # Store one energy level from the LAMDA file.
class Level:  # Level table row.
    index: int  # One-based LAMDA level index.
    energy_k: float  # Level energy in Kelvin.
    weight: float  # Statistical weight from the data file.
    label: str  # Quantum-number label as text.


@dataclass(frozen=True)  # Store one radiative transition from the LAMDA file.
class RadiativeTransition:  # Radiative transition row.
    index: int  # One-based transition index.
    upper: int  # One-based upper level index.
    lower: int  # One-based lower level index.
    aul: float  # Einstein A coefficient in s^-1.
    freq_ghz: float  # Transition frequency in GHz.
    eu_k: float  # Upper-state energy in Kelvin as listed by LAMDA.


@dataclass(frozen=True)  # Store one collision-partner block from the LAMDA file.
class CollisionPartner:  # Collision rates for a single partner.
    description: str  # Raw LAMDA partner description.
    temperatures: np.ndarray  # Temperatures where downward rates are tabulated.
    rates: dict  # Mapping (upper_level, lower_level) -> downward rates over temperatures.


@dataclass(frozen=True)  # Store all parsed information for one molecule.
class MoleculeData:  # Parsed molecule data.
    name: str  # Molecule label from the file.
    levels: dict  # Mapping level index -> Level.
    radiative: list  # List of RadiativeTransition entries.
    partners: list  # List of CollisionPartner entries.


def _is_comment(line):  # Decide whether a LAMDA line is a comment or section title.
    return line.strip().startswith("!")  # LAMDA marks headers with exclamation points.


def _next_data_index(lines, start):  # Find the next non-comment, non-empty line.
    index = start  # Start scanning from the requested row.
    while index < len(lines):  # Walk through the file until usable data is found.
        text = lines[index].strip()  # Normalize whitespace for the current row.
        if text and not _is_comment(text):  # Data rows are non-empty and not header rows.
            return index  # Return the row number of the usable data.
        index += 1  # Continue scanning after comments or blank rows.
    raise ValueError("Unexpected end of file while looking for a data row.")  # Signal malformed input.


def _find_marker(lines, start, required, forbidden=()):  # Find a header containing all required words.
    required = tuple(word.upper() for word in required)  # Compare markers case-insensitively.
    forbidden = tuple(word.upper() for word in forbidden)  # Compare excluded words case-insensitively.
    for index in range(start, len(lines)):  # Scan from the requested location.
        text = lines[index].upper()  # Normalize this candidate header.
        if all(word in text for word in required) and not any(word in text for word in forbidden):  # Match the requested header.
            return index  # Return the matching line number.
    raise ValueError(f"Could not find marker containing {required}.")  # Signal missing required section.


def _read_float_block(lines, start, count):  # Read a fixed number of floats across one or more rows.
    values = []  # Accumulate numeric values.
    index = start  # Track the current row.
    while len(values) < count and index < len(lines):  # Continue until enough values are collected.
        if not _is_comment(lines[index]):  # Skip LAMDA header rows.
            values.extend(float(value) for value in lines[index].split())  # Add all floats on this row.
        index += 1  # Move to the next row.
    if len(values) != count:  # Validate the requested number of values was read.
        raise ValueError(f"Expected {count} floats, found {len(values)}.")  # Signal malformed numeric block.
    return np.array(values, dtype=float), index  # Return values and the next unread row.


def _parse_rate_row(lines, start, n_temps):  # Read one collision-rate row, allowing wrapped rates.
    tokens = []  # Accumulate row tokens.
    index = start  # Track current file row.
    while len(tokens) < 3 + n_temps and index < len(lines):  # Need trans, upper, lower, and all rates.
        if not _is_comment(lines[index]):  # Skip header rows if encountered.
            tokens.extend(lines[index].split())  # Add tokens from this row.
        index += 1  # Advance to the next row.
    if len(tokens) < 3 + n_temps:  # Check for an incomplete row.
        raise ValueError("Incomplete collision-rate row.")  # Signal malformed collision data.
    upper = int(tokens[1])  # LAMDA column 2 is upper level.
    lower = int(tokens[2])  # LAMDA column 3 is lower level.
    rates = np.array([float(value) for value in tokens[3 : 3 + n_temps]], dtype=float)  # Rate coefficients.
    return upper, lower, rates, index  # Return parsed row fields and next unread row.


def _parse_lamda_file(path):  # Parse the subset of LAMDA needed for critical density.
    lines = path.read_text().splitlines()  # Load the ASCII data file.
    name_index = _next_data_index(lines, _find_marker(lines, 0, ("MOLECULE",)) + 1)  # Find molecule name.
    name = lines[name_index].strip()  # Store the molecule name text.
    n_level_index = _next_data_index(lines, _find_marker(lines, name_index, ("NUMBER", "ENERGY", "LEVEL")) + 1)  # Find level count.
    n_levels = int(lines[n_level_index].split()[0])  # Read the number of levels.
    level_start = _next_data_index(lines, n_level_index + 1)  # Locate the first level row.
    levels = {}  # Build a mapping from level index to Level.
    for row in lines[level_start : level_start + n_levels]:  # Read all level rows.
        parts = row.split()  # Split the level row into columns.
        level_index = int(parts[0])  # Level index is column 1.
        energy_k = float(parts[1]) * CM_TO_K  # Level energy is column 2 in cm^-1.
        weight = float(parts[2])  # Statistical weight is column 3.
        label = " ".join(parts[3:]).strip('"')  # Remaining columns describe quantum numbers.
        levels[level_index] = Level(level_index, energy_k, weight, label)  # Store this level.
    n_rad_index = _next_data_index(lines, _find_marker(lines, level_start + n_levels, ("NUMBER", "RADIATIVE", "TRANS")) + 1)  # Find radiative count.
    n_rad = int(lines[n_rad_index].split()[0])  # Read the number of radiative transitions.
    rad_start = _next_data_index(lines, n_rad_index + 1)  # Locate the first radiative row.
    radiative = []  # Build the radiative-transition list.
    for row in lines[rad_start : rad_start + n_rad]:  # Read all radiative rows.
        parts = row.split()  # Split the radiative row into columns.
        radiative.append(RadiativeTransition(int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])))  # Store needed columns.
    n_partner_marker = _find_marker(lines, rad_start + n_rad, ("NUMBER", "COLL", "PARTNER"))  # Find collision-partner count.
    n_partner_index = _next_data_index(lines, n_partner_marker + 1)  # Locate the partner-count value.
    n_partners = int(lines[n_partner_index].split()[0])  # Read the number of collision partners.
    partners = []  # Build the collision-partner list.
    cursor = n_partner_index + 1  # Continue parsing after the partner-count line.
    for _ in range(n_partners):  # Parse each collision-partner block.
        collision_marker = _find_marker(lines, cursor, ("COLLISIONS", "BETWEEN"))  # Find the partner description section.
        description_index = _next_data_index(lines, collision_marker + 1)  # Locate the description row.
        description = lines[description_index].strip()  # Store the raw partner description.
        n_coll_marker = _find_marker(lines, description_index + 1, ("NUMBER", "COLL", "TRANS"), ("TEMP", "PARTNER"))  # Find collision-row count.
        n_coll_index = _next_data_index(lines, n_coll_marker + 1)  # Locate the collision-row count value.
        n_coll = int(lines[n_coll_index].split()[0])  # Read the number of collision transitions.
        n_temp_marker = _find_marker(lines, n_coll_index + 1, ("NUMBER", "TEMP"))  # Find temperature count.
        n_temp_index = _next_data_index(lines, n_temp_marker + 1)  # Locate the temperature-count value.
        n_temp = int(lines[n_temp_index].split()[0])  # Read the number of temperatures.
        temp_marker = _find_marker(lines, n_temp_index + 1, ("TEMP",), ("NUMBER",))  # Find the temperature values section.
        temp_start = _next_data_index(lines, temp_marker + 1)  # Locate the first temperature row.
        temperatures, cursor = _read_float_block(lines, temp_start, n_temp)  # Read all temperatures.
        rate_start = _next_data_index(lines, cursor)  # Locate the first collision-rate row.
        rates = {}  # Store rates by downward transition.
        cursor = rate_start  # Begin reading rate rows.
        for _ in range(n_coll):  # Read every collision-rate row.
            upper, lower, values, cursor = _parse_rate_row(lines, cursor, n_temp)  # Parse one rate row.
            rates[(upper, lower)] = values  # Store downward rates for this pair.
        partners.append(CollisionPartner(description, temperatures, rates))  # Store this partner block.
    return MoleculeData(name, levels, radiative, partners)  # Return the parsed molecule.


def _data_path(molecule):  # Resolve a user molecule name to a data-file path.
    key = molecule.upper()  # Normalize molecule names for lookup.
    if key not in MOLECULE_FILES:  # Reject unknown molecules early.
        choices = ", ".join(sorted(MOLECULE_FILES))  # Prepare a helpful list of choices.
        raise ValueError(f"Unsupported molecule {molecule!r}; choose one of: {choices}.")  # Report valid options.
    return DATA_DIR / MOLECULE_FILES[key]  # Return the absolute path to the data file.


def load_data(molecule):  # Public helper for reading a molecule data file.
    return _parse_lamda_file(_data_path(molecule))  # Parse and return structured LAMDA data.


def _interp_rate(temperatures, values, temperature):  # Interpolate one rate coefficient to Tkin.
    if temperature < temperatures[0] or temperature > temperatures[-1]:  # Keep behavior explicit outside the table range.
        raise ValueError(f"T={temperature} K is outside the tabulated range {temperatures[0]}-{temperatures[-1]} K.")  # Avoid silent extrapolation.
    return float(np.interp(temperature, temperatures, values))  # Perform linear interpolation, exact at grid points.


def _is_electron_partner(description):  # Detect electron collision blocks.
    text = f" {description.lower()} "  # Pad text to make simple token checks safer.
    return "-e" in text or " e " in text or "electron" in text  # Exclude electron rates for Shirley H2 critical densities.


def _partner_fractions(partners, op_ratio):  # Assign fractional weights to collision partners.
    molecular_partners = [partner for partner in partners if not _is_electron_partner(partner.description)]  # Use neutral H2/He-like partners only.
    if not molecular_partners:  # Guard against data files that contain only electrons.
        raise ValueError("No neutral molecular collision partner was found.")  # Shirley table uses molecular partners.
    para = [partner for partner in molecular_partners if "ph2" in partner.description.lower() or "para" in partner.description.lower()]  # Detect para-H2.
    ortho = [partner for partner in molecular_partners if "oh2" in partner.description.lower() or "ortho" in partner.description.lower()]  # Detect ortho-H2.
    if para and ortho:  # If both H2 spin species are available, mix rates before dividing A.
        total = op_ratio + 1.0  # Convert o/p ratio into fractions.
        return [(partner, 1.0 / total) for partner in para] + [(partner, op_ratio / total) for partner in ortho]  # Return para and ortho fractions.
    fraction = 1.0 / len(molecular_partners)  # Otherwise average all neutral molecular partner blocks equally.
    return [(partner, fraction) for partner in molecular_partners]  # Return the default partner fractions.


def _downward_sum(level, temperature, partner):  # Sum all downward collisions out of the upper level.
    total = 0.0  # Initialize the downward collisional depopulation rate.
    for (upper, lower), values in partner.rates.items():  # Walk through all tabulated downward rates.
        if upper == level and lower != level:  # Select transitions from the target level to lower levels.
            total += _interp_rate(partner.temperatures, values, temperature)  # Add the interpolated downward rate.
    return total  # Return the total downward rate out of the level.


def _upward_sum(level, temperature, partner, levels):  # Sum upward collisions out of the upper level using detailed balance.
    total = 0.0  # Initialize the upward collisional depopulation rate.
    source = levels[level]  # Cache the source level object.
    for (upper, lower), values in partner.rates.items():  # Walk through all tabulated downward rates.
        if lower == level and upper != level:  # A downward upper->level row defines an upward level->upper rate.
            target = levels[upper]  # Cache the higher target level object.
            downward = _interp_rate(partner.temperatures, values, temperature)  # Interpolate gamma_upper,level.
            delta_e = target.energy_k - source.energy_k  # Energy gap in Kelvin.
            detailed_balance = target.weight / source.weight * np.exp(-delta_e / temperature)  # Shirley equation (5).
            total += downward * detailed_balance  # Add the upward rate out of the source level.
    return total  # Return the total upward rate out of the level.


def _mixed_collision_sum(level, temperature, partners, levels, op_ratio):  # Sum collision rates after partner mixing.
    total = 0.0  # Initialize the mixed collisional depopulation rate.
    for partner, fraction in _partner_fractions(partners, op_ratio):  # Apply each partner fraction to its rate sum.
        gamma_down = _downward_sum(level, temperature, partner)  # Sum downward rates for this partner.
        gamma_up = _upward_sum(level, temperature, partner, levels)  # Sum upward rates for this partner.
        total += fraction * (gamma_down + gamma_up)  # Mix rates, not already-inverted critical densities.
    return total  # Return the partner-weighted total gamma.


def _transition_from_request(data, j_low=None, transition_index=None, upper=None, lower=None):  # Resolve the target radiative transition.
    if transition_index is not None:  # Let users address complex spectra by LAMDA transition index.
        for transition in data.radiative:  # Search radiative transitions.
            if transition.index == transition_index:  # Match the requested transition index.
                return transition  # Return the requested transition.
        raise ValueError(f"Transition index {transition_index} was not found.")  # Signal invalid transition index.
    if upper is not None and lower is not None:  # Let users address complex spectra by level indices.
        for transition in data.radiative:  # Search radiative transitions.
            if transition.upper == upper and transition.lower == lower:  # Match the requested levels.
                return transition  # Return the requested transition.
        raise ValueError(f"Radiative transition upper={upper}, lower={lower} was not found.")  # Signal invalid level pair.
    if j_low is None:  # The legacy J_low argument is required when no explicit transition is supplied.
        raise ValueError("Provide J_low, transition_index, or both upper and lower.")  # Explain the required selector.
    legacy_upper = int(j_low) + 2  # Legacy API assumes LAMDA level index J+1 and transition J+1 -> J.
    legacy_lower = int(j_low) + 1  # Convert zero-based lower J into one-based lower level.
    return _transition_from_request(data, upper=legacy_upper, lower=legacy_lower)  # Resolve the legacy simple-rotor transition.


def ncrit(molecule, J_low=None, T=None, op_ratio=3.0, verbose=False, transition_index=None, upper=None, lower=None):  # Calculate Shirley nthin,nobg_crit.
    if T is None:  # Preserve a clear error when kinetic temperature is omitted.
        raise ValueError("T must be supplied in Kelvin.")  # Tkin is needed for interpolation and detailed balance.
    data = load_data(molecule)  # Parse the requested molecule file.
    transition = _transition_from_request(data, J_low, transition_index, upper, lower)  # Resolve the target radiative line.
    gamma_total = _mixed_collision_sum(transition.upper, float(T), data.partners, data.levels, float(op_ratio))  # Sum mixed collisional depopulation.
    if gamma_total <= 0.0:  # Guard against missing collision coverage for this level.
        raise ValueError(f"No collisional depopulation rates found for upper level {transition.upper}.")  # Report unusable transition.
    value = transition.aul / gamma_total  # Shirley equation (4): ncrit = A_jk / sum_i gamma_ji.
    if verbose:  # Optionally print a compact report.
        print(f"{molecule} transition {transition.upper}->{transition.lower}")  # Identify the transition by LAMDA levels.
        print(f"Tkin = {float(T):.3g} K")  # Report kinetic temperature.
        print(f"ortho/para H2 ratio = {float(op_ratio):.3g}")  # Report the requested o/p ratio.
        print(f"Aul = {transition.aul:.6e} s^-1")  # Report Einstein A.
        print(f"sum gamma = {gamma_total:.6e} cm^3 s^-1")  # Report the denominator.
        print(f"critical density = {value:.6e} cm^-3")  # Report the final critical density.
    return value  # Return the scalar critical density.


def input(molecule):  # Legacy helper kept for old notebooks that imported input().
    data = load_data(molecule)  # Parse the molecule file.
    neutral_partners = [partner for partner in data.partners if not _is_electron_partner(partner.description)]  # Ignore electron blocks.
    para_partners = [partner for partner in neutral_partners if "ph2" in partner.description.lower() or "para" in partner.description.lower()]  # Find para-H2 blocks.
    ortho_partners = [partner for partner in neutral_partners if "oh2" in partner.description.lower() or "ortho" in partner.description.lower()]  # Find ortho-H2 blocks.
    para_partner = para_partners[0] if para_partners else neutral_partners[0]  # Use para-H2 when present, otherwise the first neutral partner.
    ortho_partner = ortho_partners[0] if ortho_partners else para_partner  # Use ortho-H2 when present, otherwise mirror para.
    temperatures = para_partner.temperatures  # Return the para/neutral temperature grid.
    para_collisions = _legacy_collision_array(para_partner)  # Convert para/neutral rates to the old array layout.
    ortho_collisions = _legacy_collision_array(ortho_partner)  # Convert ortho/neutral rates to the old array layout.
    a_values = np.array([transition.aul for transition in data.radiative], dtype=float)  # Return Einstein A values.
    freq_values = np.array([transition.freq_ghz * 1e9 for transition in data.radiative], dtype=float)  # Return Hz frequencies.
    eu_values = np.array([0.0] + [data.levels[index].energy_k for index in sorted(data.levels)], dtype=float)  # Return old one-offset level energies.
    return temperatures, para_collisions, ortho_collisions, a_values, freq_values, eu_values  # Match the old tuple shape.


def _legacy_collision_array(partner):  # Convert parsed collision rates to the old TRANS, UP, LOW, RATE... layout.
    rows = []  # Accumulate legacy rows.
    for row_index, ((upper, lower), values) in enumerate(partner.rates.items(), start=1):  # Walk stored rates.
        rows.append([row_index, upper, lower, *values])  # Rebuild one legacy collision row.
    return np.array(rows, dtype=float)  # Return the legacy collision array.
