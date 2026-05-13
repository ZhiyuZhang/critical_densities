"""NH3 HFS critical-density demos using Loreau et al. data files."""

import os
from pathlib import Path

import numpy as np

from critical_density_hfs import compute_hfs_ncrit, load_hfs_molecule


BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "outputs"
DEMO_TK = 50.0
DEMO_PARTNER = "p-H2"
LTE_TEX = 10.0
GAUSSIAN_FWHM_KMS = 0.20
C_KMS = 299792.458
P_NH3_DATA = BASE_DIR / "p-nh3@loreau.dat"
O_NH3_DATA = BASE_DIR / "o-nh3@loreau.dat"


def _value_label(value):
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _level_label(level):
    qn = level.qn
    if {"J", "K", "F"}.issubset(qn):
        parity_key = "sym" if "sym" in qn else "eps" if "eps" in qn else None
        parity = f",{parity_key}={_value_label(qn[parity_key])}" if parity_key else ""
        return f"J={_value_label(qn['J'])},K={_value_label(qn['K'])}{parity},F={_value_label(qn['F'])}"
    return level.label


def _compute_rows(molecule):
    level_rows, _ = compute_hfs_ncrit(
        molecule,
        [DEMO_TK],
        partner=DEMO_PARTNER,
        require_hfs_collisions=True,
    )
    return {row["upper_id"]: row for row in level_rows}


def _select_p_nh3_11_inversion(molecule):
    # p-NH3 (1,1) inversion components connect J=1,K=1,sym=-1 upper levels to sym=+1 lower levels.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 1
        and molecule.levels[transition.upper_id].qn.get("K") == 1
        and molecule.levels[transition.upper_id].qn.get("sym") == -1
        and molecule.levels[transition.lower_id].qn.get("J") == 1
        and molecule.levels[transition.lower_id].qn.get("K") == 1
        and molecule.levels[transition.lower_id].qn.get("sym") == 1
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _select_p_nh3_22_inversion(molecule):
    # p-NH3 (2,2) inversion components connect J=2,K=2,sym=-1 upper levels to sym=+1 lower levels.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 2
        and molecule.levels[transition.upper_id].qn.get("K") == 2
        and molecule.levels[transition.upper_id].qn.get("sym") == -1
        and molecule.levels[transition.lower_id].qn.get("J") == 2
        and molecule.levels[transition.lower_id].qn.get("K") == 2
        and molecule.levels[transition.lower_id].qn.get("sym") == 1
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _select_p_nh3_21_11_rotational(molecule):
    # p-NH3 2_1-1_1 rotational-inversion branches appear as two parity branches in this file.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 2
        and molecule.levels[transition.upper_id].qn.get("K") == 1
        and molecule.levels[transition.lower_id].qn.get("J") == 1
        and molecule.levels[transition.lower_id].qn.get("K") == 1
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _select_p_nh3_22_11_if_present(molecule):
    # This checks the user's possible 2_2-1_1 notation; the Loreau file has no such radiative branch.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 2
        and molecule.levels[transition.upper_id].qn.get("K") == 2
        and molecule.levels[transition.lower_id].qn.get("J") == 1
        and molecule.levels[transition.lower_id].qn.get("K") == 1
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _select_o_nh3_10_00(molecule):
    # o-NH3 lowest rotational HFS components connect J=1,K=0,eps=+1 upper levels to J=0,K=0,eps=-1.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 1
        and molecule.levels[transition.upper_id].qn.get("K") == 0
        and molecule.levels[transition.upper_id].qn.get("eps") == 1
        and molecule.levels[transition.lower_id].qn.get("J") == 0
        and molecule.levels[transition.lower_id].qn.get("K") == 0
        and molecule.levels[transition.lower_id].qn.get("eps") == -1
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _print_table(title, transitions, molecule, ncrit_by_upper):
    if not transitions:
        print(f"\n{title}: no matching radiative transitions in this data file.")
        return
    print(f"\n{title} at Tk={DEMO_TK:g} K, partner={DEMO_PARTNER}")
    print("freq_GHz    upper_id  lower_id  upper_level                         lower_level                         ncrit_cm-3")
    for transition in transitions:
        upper = molecule.levels[transition.upper_id]
        lower = molecule.levels[transition.lower_id]
        ncrit = ncrit_by_upper[transition.upper_id]["ncrit_cm-3"]
        print(
            f"{transition.frequency_GHz:10.5f}  "
            f"{transition.upper_id:8d}  "
            f"{transition.lower_id:8d}  "
            f"{_level_label(upper):35s} "
            f"{_level_label(lower):35s} "
            f"{ncrit:10.3e}"
        )


def _plot(title, filename, transitions, molecule, ncrit_by_upper):
    if not transitions:
        return
    cache_dir = BASE_DIR / ".plot_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    labels = []
    values = []
    for transition in transitions:
        upper = molecule.levels[transition.upper_id]
        lower = molecule.levels[transition.lower_id]
        labels.append(f"{transition.frequency_GHz:.3f} GHz\nu:{_level_label(upper)}\nl:{_level_label(lower)}")
        values.append(ncrit_by_upper[transition.upper_id]["ncrit_cm-3"])

    fig_width = max(12, 1.25 * len(values))
    fig, ax = plt.subplots(figsize=(fig_width, 7), constrained_layout=True)
    x = np.arange(len(values))
    bars = ax.bar(x, values, color="#E15759", edgecolor="#6F2526", linewidth=0.7)
    ax.set_yscale("log")
    ax.set_ylabel(r"$n_{\rm crit}$ (cm$^{-3}$)")
    ax.set_xlabel("NH3 hyperfine transition, sorted by frequency")
    ax.set_title(f"{title} level-based critical densities at Tk={DEMO_TK:g} K ({DEMO_PARTNER})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=7)
    ax.grid(axis="y", which="both", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.05, f"{value:.1e}", ha="center", va="bottom", fontsize=7, rotation=90)
    out_path = OUTPUT_DIR / filename
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved NH3 ncrit histogram: {out_path}")


def _lte_component_weight(transition, molecule):
    upper = molecule.levels[transition.upper_id]
    # Optically thin LTE integrated intensity is proportional to upper population times A_ul.
    return upper.g * np.exp(-upper.energy_K / LTE_TEX) * transition.A_s


def _velocity_offset_kms(frequency_GHz, reference_GHz):
    return C_KMS * (reference_GHz - frequency_GHz) / reference_GHz


def _plot_lte_spectrum(title, filename, transitions, molecule, ncrit_by_upper):
    if not transitions:
        return
    cache_dir = BASE_DIR / ".plot_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    reference = np.mean([transition.frequency_GHz for transition in transitions])
    centers = np.array([_velocity_offset_kms(transition.frequency_GHz, reference) for transition in transitions])
    weights = np.array([_lte_component_weight(transition, molecule) for transition in transitions])
    weights = weights / weights.max()
    sigma = GAUSSIAN_FWHM_KMS / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    x_min = centers.min() - 5.0 * GAUSSIAN_FWHM_KMS
    x_max = centers.max() + 5.0 * GAUSSIAN_FWHM_KMS
    velocity = np.linspace(x_min, x_max, 4000)

    upper_ids = sorted({transition.upper_id for transition in transitions})
    colors = dict(zip(upper_ids, plt.cm.tab10(np.linspace(0, 1, max(len(upper_ids), 2)))[: len(upper_ids)]))
    component_profiles = []
    total = np.zeros_like(velocity)
    for transition, center, weight in zip(transitions, centers, weights):
        profile = weight * np.exp(-0.5 * ((velocity - center) / sigma) ** 2)
        component_profiles.append(profile)
        total += profile
    norm = total.max()
    total = total / norm
    component_profiles = [profile / norm for profile in component_profiles]

    fig, ax = plt.subplots(figsize=(14, 7), constrained_layout=True)
    for transition, center, profile in zip(transitions, centers, component_profiles):
        color = colors[transition.upper_id]
        ax.fill_between(velocity, profile, color=color, alpha=0.28)
        ax.plot(velocity, profile, color=color, lw=1.0, alpha=0.85)
        ax.axvline(center, color=color, lw=0.8, alpha=0.65)
    ax.plot(velocity, total, color="black", lw=2.0, label="LTE optically thin sum")

    handles = [plt.Line2D([0], [0], color="black", lw=2.0, label="LTE optically thin sum")]
    for upper_id in upper_ids:
        level = molecule.levels[upper_id]
        ncrit = ncrit_by_upper[upper_id]["ncrit_cm-3"]
        handles.append(
            plt.Line2D(
                [0],
                [0],
                color=colors[upper_id],
                lw=5.0,
                label=f"u={_level_label(level)}; ncrit={ncrit:.2e} cm^-3",
            )
        )
    ax.legend(handles=handles, fontsize=8, loc="upper right", frameon=False)
    ax.set_xlabel(f"Velocity offset relative to {reference:.6f} GHz (km s$^{{-1}}$)")
    ax.set_ylabel("Normalized intensity")
    ax.set_title(
        f"{title}: Gaussian LTE optically thin HFS spectrum\n"
        f"Tex={LTE_TEX:g} K, FWHM={GAUSSIAN_FWHM_KMS:g} km s$^{{-1}}$, ncrit at Tk={DEMO_TK:g} K"
    )
    ax.grid(alpha=0.25)
    ax.invert_xaxis()
    out_path = OUTPUT_DIR / filename
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved NH3 LTE HFS spectrum: {out_path}")


def main():
    p_molecule = load_hfs_molecule(P_NH3_DATA, "p-NH3")
    p_ncrit = _compute_rows(p_molecule)
    p11_transitions = _select_p_nh3_11_inversion(p_molecule)
    _print_table("p-NH3 (1,1) inversion HFS components", p11_transitions, p_molecule, p_ncrit)
    _plot("p-NH3 (1,1) inversion HFS", "p_nh3_11_inversion_ncrit_hist.png", p11_transitions, p_molecule, p_ncrit)
    _plot_lte_spectrum("p-NH3 (1,1) inversion HFS", "p_nh3_11_inversion_lte_spectrum_ncrit.png", p11_transitions, p_molecule, p_ncrit)

    p22_transitions = _select_p_nh3_22_inversion(p_molecule)
    _print_table("p-NH3 (2,2) inversion HFS components", p22_transitions, p_molecule, p_ncrit)
    _plot("p-NH3 (2,2) inversion HFS", "p_nh3_22_inversion_ncrit_hist.png", p22_transitions, p_molecule, p_ncrit)
    _plot_lte_spectrum("p-NH3 (2,2) inversion HFS", "p_nh3_22_inversion_lte_spectrum_ncrit.png", p22_transitions, p_molecule, p_ncrit)

    p21_11_transitions = _select_p_nh3_21_11_rotational(p_molecule)
    _print_table("p-NH3 2_1-1_1 rotational-inversion HFS components", p21_11_transitions, p_molecule, p_ncrit)
    _plot("p-NH3 2_1-1_1 rotational-inversion HFS", "p_nh3_21_11_ncrit_hist.png", p21_11_transitions, p_molecule, p_ncrit)
    _plot_lte_spectrum("p-NH3 2_1-1_1 rotational-inversion HFS", "p_nh3_21_11_lte_spectrum_ncrit.png", p21_11_transitions, p_molecule, p_ncrit)

    p22_11_transitions = _select_p_nh3_22_11_if_present(p_molecule)
    _print_table("p-NH3 2_2-1_1 HFS components", p22_11_transitions, p_molecule, p_ncrit)

    o_molecule = load_hfs_molecule(O_NH3_DATA, "o-NH3")
    o_transitions = _select_o_nh3_10_00(o_molecule)
    o_ncrit = _compute_rows(o_molecule)
    _print_table("o-NH3 1_0-0_0 rotational HFS components", o_transitions, o_molecule, o_ncrit)
    _plot("o-NH3 1_0-0_0 rotational HFS", "o_nh3_10_00_ncrit_hist.png", o_transitions, o_molecule, o_ncrit)


if __name__ == "__main__":
    main()
