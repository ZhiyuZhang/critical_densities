"""HFS-resolved optically thin critical-density calculator."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


CM_TO_K = 1.438776877


@dataclass(frozen=True)
class HFSLevel:
    id: int
    energy_K: float
    g: float
    label: str
    qn: dict
    is_hfs_resolved: bool


@dataclass(frozen=True)
class RadiativeTransition:
    index: int
    upper_id: int
    lower_id: int
    A_s: float
    frequency_GHz: float
    E_u_K: float


@dataclass(frozen=True)
class CollisionPartner:
    name: str
    description: str
    temperatures_K: np.ndarray
    rates: dict
    is_hfs_resolved: bool
    is_approximate: bool


@dataclass(frozen=True)
class HFSMolecule:
    species: str
    datafile: Path
    levels: dict
    radiative: list
    partners: list
    data_mode: str
    is_hfs_resolved: bool


def _is_comment(line: str) -> bool:
    return line.strip().startswith("!")


def _next_data_index(lines: list[str], start: int) -> int:
    for index in range(start, len(lines)):
        text = lines[index].strip()
        if text and not _is_comment(text):
            return index
    raise ValueError("Unexpected end of file while looking for a data row.")


def _find_marker(lines: list[str], start: int, required: tuple[str, ...], forbidden: tuple[str, ...] = ()) -> int:
    required_upper = tuple(word.upper() for word in required)
    forbidden_upper = tuple(word.upper() for word in forbidden)
    for index in range(start, len(lines)):
        text = lines[index].upper()
        if all(word in text for word in required_upper) and not any(word in text for word in forbidden_upper):
            return index
    raise ValueError(f"Could not find marker containing {required}.")


def _read_float_block(lines: list[str], start: int, count: int) -> tuple[np.ndarray, int]:
    values = []
    index = start
    while len(values) < count and index < len(lines):
        if not _is_comment(lines[index]):
            values.extend(float(value) for value in lines[index].split())
        index += 1
    if len(values) != count:
        raise ValueError(f"Expected {count} floats, found {len(values)}.")
    return np.array(values, dtype=float), index


def _parse_qn_names(level_header: str) -> list[str]:
    if "+" in level_header:
        right_side = level_header.split("+")[-1]
    elif "|" in level_header:
        right_side = level_header.split("|")[-1]
    else:
        return []
    right_side = right_side.replace("(", " ").replace(")", " ")
    names = [token for token in right_side.split() if token.upper() not in {"QNUM", "QNS"}]
    return names


def _parse_qn_dict(qn_names: list[str], qn_values: list[str]) -> dict:
    qn = {}
    if len(qn_names) == 1 and "_" in qn_names[0] and len(qn_values) == 1:
        names = qn_names[0].split("_")
        values = qn_values[0].split("_")
        return {name: _format_qn_value(value) for name, value in zip(names, values)}
    if len(qn_names) == 1 and "_" not in qn_names[0] and len(qn_values) == 1 and "_" in qn_values[0]:
        return {qn_names[0]: qn_values[0]}
    for name, value in zip(qn_names, qn_values):
        if "_" in name and "_" in value:
            for sub_name, sub_value in zip(name.split("_"), value.split("_")):
                qn[sub_name] = _format_qn_value(sub_value)
        else:
            qn[name] = _format_qn_value(value)
    return qn


def _format_qn_value(value: str):
    try:
        number = float(value)
    except ValueError:
        return value.strip('"')
    return int(number) if number.is_integer() else number


def _partner_name(description: str) -> str:
    lowered = description.lower()
    if "para" in lowered or "p-h2" in lowered or "ph2" in lowered:
        return "p-H2"
    if "ortho" in lowered or "o-h2" in lowered or "oh2" in lowered:
        return "o-H2"
    if "h2" in lowered:
        return "H2"
    if "electron" in lowered or "-e" in lowered:
        return "e"
    if "he" in lowered:
        return "He"
    return description.split(maxsplit=1)[-1]


def _format_qn_number(value) -> str:
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _level_group(level: HFSLevel) -> str:
    # Group HFS sublevels by the parent rotational/fine-structure quantum number for easier scanning.
    if "N" in level.qn:
        return f"N={_format_qn_number(level.qn['N'])}"
    if "J" in level.qn and "K" in level.qn and ("sym" in level.qn or "eps" in level.qn):
        parity_key = "sym" if "sym" in level.qn else "eps"
        return f"J={_format_qn_number(level.qn['J'])},K={_format_qn_number(level.qn['K'])},{parity_key}={_format_qn_number(level.qn[parity_key])}"
    if "J" in level.qn and "K" in level.qn:
        return f"J={_format_qn_number(level.qn['J'])},K={_format_qn_number(level.qn['K'])}"
    if "J" in level.qn:
        return f"J={_format_qn_number(level.qn['J'])}"
    return f"level={level.id}"


def _parse_rate_row(lines: list[str], start: int, n_temps: int) -> tuple[int, int, np.ndarray, int]:
    tokens = []
    index = start
    while len(tokens) < 3 + n_temps and index < len(lines):
        if not _is_comment(lines[index]):
            tokens.extend(lines[index].split())
        index += 1
    if len(tokens) < 3 + n_temps:
        raise ValueError("Incomplete collision-rate row.")
    upper = int(tokens[1])
    lower = int(tokens[2])
    rates = np.array([float(value) for value in tokens[3 : 3 + n_temps]], dtype=float)
    return upper, lower, rates, index


def load_hfs_molecule(datafile: str | Path, species: str | None = None) -> HFSMolecule:
    path = Path(datafile).expanduser().resolve()
    lines = path.read_text().splitlines()
    name_index = _next_data_index(lines, _find_marker(lines, 0, ("MOLECULE",)) + 1)
    molecule_name = species or lines[name_index].strip()
    n_level_index = _next_data_index(lines, _find_marker(lines, name_index, ("NUMBER", "ENERGY", "LEVEL")) + 1)
    n_levels = int(lines[n_level_index].split()[0])
    level_header_index = _find_marker(lines, n_level_index + 1, ("LEVEL",))
    qn_names = _parse_qn_names(lines[level_header_index])
    level_start = _next_data_index(lines, level_header_index + 1)
    levels = {}
    for row in lines[level_start : level_start + n_levels]:
        parts = row.split()
        level_id = int(parts[0])
        qn_values = parts[3:]
        qn = _parse_qn_dict(qn_names, qn_values)
        label = ",".join(f"{key}={value}" for key, value in qn.items()) or " ".join(qn_values).strip('"')
        # Presence of F in the level quantum numbers is the most robust marker for HFS-resolved levels.
        is_hfs = any(key.upper() == "F" for key in qn)
        levels[level_id] = HFSLevel(level_id, float(parts[1]) * CM_TO_K, float(parts[2]), label, qn, is_hfs)

    n_rad_index = _next_data_index(lines, _find_marker(lines, level_start + n_levels, ("NUMBER", "RADIATIVE", "TRANS")) + 1)
    n_rad = int(lines[n_rad_index].split()[0])
    rad_start = _next_data_index(lines, n_rad_index + 1)
    radiative = []
    for row in lines[rad_start : rad_start + n_rad]:
        parts = row.split()
        radiative.append(RadiativeTransition(int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])))

    n_partner_index = _next_data_index(lines, _find_marker(lines, rad_start + n_rad, ("NUMBER", "COLL", "PARTNER")) + 1)
    n_partners = int(lines[n_partner_index].split()[0])
    partners = []
    cursor = n_partner_index + 1
    for _ in range(n_partners):
        try:
            marker = _find_marker(lines, cursor, ("COLLISIONS", "BETWEEN"))
        except ValueError:
            marker = _find_marker(lines, cursor, ("COLLISIONS", "WITH"))
        desc_index = _next_data_index(lines, marker + 1)
        description = lines[desc_index].strip()
        n_coll_index = _next_data_index(lines, _find_marker(lines, desc_index + 1, ("NUMBER", "COLL", "TRANS"), ("TEMP",)) + 1)
        n_coll = int(lines[n_coll_index].split()[0])
        n_temp_index = _next_data_index(lines, _find_marker(lines, n_coll_index + 1, ("NUMBER", "TEMP")) + 1)
        n_temp = int(lines[n_temp_index].split()[0])
        temp_start = _next_data_index(lines, _find_marker(lines, n_temp_index + 1, ("TEMP",), ("NUMBER",)) + 1)
        temperatures, cursor = _read_float_block(lines, temp_start, n_temp)
        rate_start = _next_data_index(lines, cursor)
        rates = {}
        cursor = rate_start
        for _ in range(n_coll):
            upper, lower, values, cursor = _parse_rate_row(lines, cursor, n_temp)
            rates[(upper, lower)] = values
        all_ids_are_hfs = all(levels[upper].is_hfs_resolved and levels[lower].is_hfs_resolved for upper, lower in rates)
        is_approximate = "approx" in description.lower() or "random" in description.lower()
        partners.append(CollisionPartner(_partner_name(description), description, temperatures, rates, all_ids_are_hfs, is_approximate))

    is_hfs_resolved = any(level.is_hfs_resolved for level in levels.values())
    collision_hfs = any(partner.is_hfs_resolved for partner in partners)
    if is_hfs_resolved and collision_hfs:
        data_mode = "HFS-resolved spectroscopy + HFS-resolved collisions"
    elif is_hfs_resolved:
        data_mode = "HFS-resolved spectroscopy only"
    else:
        data_mode = "non-HFS / collapsed"
    return HFSMolecule(molecule_name, path, levels, radiative, partners, data_mode, is_hfs_resolved)


def _interp_rate(temperatures: np.ndarray, values: np.ndarray, temperature: float, allow_extrapolation: bool) -> tuple[float, str | None]:
    warning = None
    if temperature < temperatures[0] or temperature > temperatures[-1]:
        if not allow_extrapolation:
            raise ValueError(f"T={temperature} K is outside the tabulated range {temperatures[0]}-{temperatures[-1]} K.")
        warning = f"extrapolated from {temperatures[0]}-{temperatures[-1]} K"

    # Collision rates vary approximately exponentially with log(T); use log-log only when rates are strictly positive.
    if np.all(values > 0):
        x = np.log(temperatures)
        y = np.log(values)
        target = np.log(temperature)
        result = _linear_interp_or_extrap(x, y, target)
        return float(np.exp(result)), warning
    result = _linear_interp_or_extrap(temperatures, values, temperature)
    return max(float(result), 0.0), warning


def _linear_interp_or_extrap(x: np.ndarray, y: np.ndarray, target: float) -> float:
    if x[0] <= target <= x[-1]:
        return float(np.interp(target, x, y))
    if target < x[0]:
        left, right = 0, 1
    else:
        left, right = -2, -1
    slope = (y[right] - y[left]) / (x[right] - x[left])
    return float(y[left] + slope * (target - x[left]))


def _select_partners(molecule: HFSMolecule, partner: str | None, h2_opr: float | None) -> list[tuple[CollisionPartner, float, str]]:
    available = {item.name: item for item in molecule.partners}
    if h2_opr is not None:
        if "p-H2" not in available or "o-H2" not in available:
            raise ValueError("h2_opr requires both p-H2 and o-H2 collision partners.")
        total = h2_opr + 1.0
        return [(available["p-H2"], 1.0 / total, f"p-H2 fraction={1.0 / total:.6g}"), (available["o-H2"], h2_opr / total, f"o-H2 fraction={h2_opr / total:.6g}")]
    if partner is None:
        raise ValueError(f"Choose a collision partner explicitly. Available: {', '.join(sorted(available))}.")
    if partner not in available:
        raise ValueError(f"Partner {partner!r} not found. Available: {', '.join(sorted(available))}.")
    return [(available[partner], 1.0, None)]


def _collision_sum_for_level(
    molecule: HFSMolecule,
    level_id: int,
    temperature: float,
    selected_partners: list[tuple[CollisionPartner, float, str | None]],
    allow_extrapolation: bool,
) -> tuple[float, int, list[str]]:
    total = 0.0
    branch_count = 0
    warnings = []
    source = molecule.levels[level_id]
    for partner, fraction, partner_note in selected_partners:
        if partner_note:
            warnings.append(partner_note)
        for (upper, lower), values in partner.rates.items():
            if upper == level_id and lower != level_id:
                rate, warning = _interp_rate(partner.temperatures_K, values, temperature, allow_extrapolation)
                total += fraction * rate
                branch_count += 1
                if warning:
                    warnings.append(warning)
            elif lower == level_id and upper != level_id:
                # LAMDA lists downward de-excitation; detailed balance reconstructs upward excitation out of this level.
                downward, warning = _interp_rate(partner.temperatures_K, values, temperature, allow_extrapolation)
                target = molecule.levels[upper]
                delta_E = target.energy_K - source.energy_K
                rate = downward * (target.g / source.g) * np.exp(-delta_E / temperature)
                total += fraction * rate
                branch_count += 1
                if warning:
                    warnings.append(warning)
    return total, branch_count, sorted(set(warnings))


def compute_hfs_ncrit(
    molecule: HFSMolecule,
    temperatures: list[float],
    partner: str | None = None,
    h2_opr: float | None = None,
    require_hfs_collisions: bool = False,
    allow_approximate_hfs_collisions: bool = False,
    allow_extrapolation: bool = False,
) -> tuple[list[dict], list[dict]]:
    selected_partners = _select_partners(molecule, partner, h2_opr)
    if require_hfs_collisions and not all(item.is_hfs_resolved for item, _, _ in selected_partners):
        raise ValueError("HFS-resolved collisions are required, but the selected collision data are not HFS-resolved.")
    if any(item.is_approximate for item, _, _ in selected_partners) and not allow_approximate_hfs_collisions:
        raise ValueError("Approximate HFS collision data require allow_approximate_hfs_collisions=True.")

    radiative_by_upper = {}
    for transition in molecule.radiative:
        if molecule.levels[transition.lower_id].energy_K < molecule.levels[transition.upper_id].energy_K:
            radiative_by_upper.setdefault(transition.upper_id, []).append(transition)

    partner_label = f"H2_OPR={h2_opr:g}" if h2_opr is not None else partner
    level_rows = []
    branch_rows = []
    for temperature in temperatures:
        ncrit_by_upper = {}
        for upper_id, branches in sorted(radiative_by_upper.items()):
            A_tot = sum(branch.A_s for branch in branches)
            gamma_tot, n_coll, warnings = _collision_sum_for_level(molecule, upper_id, float(temperature), selected_partners, allow_extrapolation)
            if gamma_tot <= 0.0:
                warnings.append("no collisional depopulation rates found")
                ncrit = np.nan
            else:
                ncrit = A_tot / gamma_tot
            upper = molecule.levels[upper_id]
            ncrit_by_upper[upper_id] = ncrit
            level_rows.append(
                {
                    "species": molecule.species,
                    "datafile": molecule.datafile.name,
                    "partner": partner_label,
                    "Tk_K": float(temperature),
                    "upper_id": upper_id,
                    "upper_group": _level_group(upper),
                    "upper_label": upper.label,
                    "upper_qn_json": json.dumps(upper.qn, sort_keys=True),
                    "E_u_K": upper.energy_K,
                    "g_u": upper.g,
                    "A_tot_s-1": A_tot,
                    "gamma_tot_cm3_s-1": gamma_tot,
                    "ncrit_cm-3": ncrit,
                    "n_rad_branches": len(branches),
                    "n_coll_branches": n_coll,
                    "data_mode": molecule.data_mode,
                    "is_hfs_resolved": molecule.is_hfs_resolved,
                    "is_collision_hfs_resolved": all(item.is_hfs_resolved for item, _, _ in selected_partners),
                    "is_collision_approximate": any(item.is_approximate for item, _, _ in selected_partners),
                    "warnings": "; ".join(warnings),
                }
            )
        for transition in molecule.radiative:
            upper_id = transition.upper_id
            if upper_id not in ncrit_by_upper:
                continue
            A_tot = sum(branch.A_s for branch in radiative_by_upper[upper_id])
            upper = molecule.levels[upper_id]
            lower = molecule.levels[transition.lower_id]
            branch_rows.append(
                {
                    "species": molecule.species,
                    "upper_id": upper_id,
                    "lower_id": transition.lower_id,
                    "upper_group": _level_group(upper),
                    "lower_group": _level_group(lower),
                    "upper_label": upper.label,
                    "lower_label": lower.label,
                    "frequency_GHz": transition.frequency_GHz,
                    "A_ul_s-1": transition.A_s,
                    "A_tot_upper_s-1": A_tot,
                    "branching_ratio": transition.A_s / A_tot if A_tot > 0 else np.nan,
                    "upper_ncrit_cm-3": ncrit_by_upper[upper_id],
                    "Tk_K": float(temperature),
                    "partner": partner_label,
                    "data_mode": molecule.data_mode,
                }
            )
    return level_rows, branch_rows


def write_csv(rows: list[dict], path: str | Path) -> None:
    path = Path(path)
    if not rows:
        raise ValueError(f"No rows to write for {path}.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def _parse_temperatures(text: str | None, grid: str | None) -> list[float]:
    if text:
        return [float(item) for item in text.split(",") if item.strip()]
    if grid:
        start, stop, count = grid.split(":")
        return [float(item) for item in np.linspace(float(start), float(stop), int(count))]
    raise ValueError("Provide --temperatures or --temperature-grid.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compute level-based HFS-resolved critical densities.")
    parser.add_argument("--data", required=True, help="LAMDA-style HFS data file.")
    parser.add_argument("--species", default=None, help="Species label to write in outputs.")
    parser.add_argument("--partner", default=None, help="Collision partner, for example p-H2 or o-H2.")
    parser.add_argument("--h2-opr", type=float, default=None, help="Use an ortho/para H2 mixture at this OPR.")
    parser.add_argument("--temperatures", default=None, help="Comma-separated temperatures, e.g. 10,20,50.")
    parser.add_argument("--temperature-grid", default=None, help="Grid Tmin:Tmax:N.")
    parser.add_argument("--require-hfs-collisions", action="store_true", help="Fail unless selected collision data are HFS-resolved.")
    parser.add_argument("--allow-approximate-hfs-collisions", action="store_true", help="Allow collision data marked approximate.")
    parser.add_argument("--allow-extrapolation", action="store_true", help="Allow temperature extrapolation and mark warnings.")
    parser.add_argument("--out", required=True, help="Output directory.")
    parser.add_argument("--verbose", action="store_true", help="Print a short summary.")
    args = parser.parse_args(argv)

    molecule = load_hfs_molecule(args.data, args.species)
    temperatures = _parse_temperatures(args.temperatures, args.temperature_grid)
    level_rows, branch_rows = compute_hfs_ncrit(
        molecule,
        temperatures,
        partner=args.partner,
        h2_opr=args.h2_opr,
        require_hfs_collisions=args.require_hfs_collisions,
        allow_approximate_hfs_collisions=args.allow_approximate_hfs_collisions,
        allow_extrapolation=args.allow_extrapolation,
    )
    out_dir = Path(args.out)
    write_csv(level_rows, out_dir / "level_ncrit.csv")
    write_csv(branch_rows, out_dir / "transition_branches.csv")
    if args.verbose:
        partners = ", ".join(f"{partner.name}: {partner.description}" for partner in molecule.partners)
        print(f"species: {molecule.species}")
        print(f"data mode: {molecule.data_mode}")
        print(f"partners: {partners}")
        print(f"wrote {len(level_rows)} level rows and {len(branch_rows)} transition rows to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
