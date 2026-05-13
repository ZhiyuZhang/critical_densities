"""Small toy checks for the HFS critical-density implementation."""

import math
import os
from pathlib import Path

import numpy as np

from critical_density_hfs import (
    CollisionPartner,
    HFSMolecule,
    HFSLevel,
    RadiativeTransition,
    _collision_sum_for_level,
    compute_hfs_ncrit,
    load_hfs_molecule,
)


BASE_DIR = Path(__file__).resolve().parent
CN_DATA = BASE_DIR / "cn-hfs.dat"
OUTPUT_DIR = BASE_DIR / "outputs"
CN_DEMO_TK = 10.0
CN_DEMO_PARTNER = "p-H2"
LTE_TEX = 10.0
GAUSSIAN_FWHM_KMS = 2.0
C_KMS = 299792.458


def toy_branching_molecule():
    levels = {
        1: HFSLevel(1, 0.0, 2.0, "N=0,J=0.5,F=0.5", {"N": 0, "J": 0.5, "F": 0.5}, True),
        2: HFSLevel(2, 1.0, 4.0, "N=0,J=0.5,F=1.5", {"N": 0, "J": 0.5, "F": 1.5}, True),
        3: HFSLevel(3, 2.0, 6.0, "N=0,J=1.5,F=2.5", {"N": 0, "J": 1.5, "F": 2.5}, True),
        4: HFSLevel(4, 10.0, 8.0, "N=1,J=1.5,F=3.5", {"N": 1, "J": 1.5, "F": 3.5}, True),
    }
    radiative = [
        RadiativeTransition(1, 4, 1, 1e-5, 100.0, 10.0),
        RadiativeTransition(2, 4, 2, 2e-5, 101.0, 10.0),
        RadiativeTransition(3, 4, 3, 3e-5, 102.0, 10.0),
    ]
    partner = CollisionPartner(
        "p-H2",
        "toy p-H2 exact HFS",
        np.array([10.0, 20.0]),
        {
            (4, 1): np.array([1e-10, 1e-10]),
            (4, 2): np.array([2e-10, 2e-10]),
            (4, 3): np.array([3e-10, 3e-10]),
        },
        True,
        False,
    )
    return HFSMolecule("toy", Path("toy.dat"), levels, radiative, [partner], "toy HFS", True)


def test_branching_and_collision_sums():
    molecule = toy_branching_molecule()
    level_rows, branch_rows = compute_hfs_ncrit(molecule, [10.0], partner="p-H2", require_hfs_collisions=True)
    assert len(level_rows) == 1
    assert math.isclose(level_rows[0]["A_tot_s-1"], 6e-5)
    assert math.isclose(level_rows[0]["gamma_tot_cm3_s-1"], 6e-10)
    assert math.isclose(level_rows[0]["ncrit_cm-3"], 1e5)
    assert {row["upper_ncrit_cm-3"] for row in branch_rows} == {level_rows[0]["ncrit_cm-3"]}
    assert math.isclose(sum(row["branching_ratio"] for row in branch_rows), 1.0)


def test_detailed_balance_reconstructs_upward_rate():
    levels = {
        1: HFSLevel(1, 0.0, 2.0, "N=0,J=0.5,F=0.5", {"N": 0, "J": 0.5, "F": 0.5}, True),
        2: HFSLevel(2, 10.0, 6.0, "N=1,J=1.5,F=2.5", {"N": 1, "J": 1.5, "F": 2.5}, True),
    }
    partner = CollisionPartner("p-H2", "toy p-H2 exact HFS", np.array([10.0]), {(2, 1): np.array([1e-10])}, True, False)
    molecule = HFSMolecule("toy", Path("toy.dat"), levels, [], [partner], "toy HFS", True)
    gamma, branches, warnings = _collision_sum_for_level(molecule, 1, 10.0, [(partner, 1.0, None)], False)
    expected = 1e-10 * (6.0 / 2.0) * math.exp(-10.0 / 10.0)
    assert branches == 1
    assert warnings == []
    assert math.isclose(gamma, expected)


def _half_label(value):
    if isinstance(value, float) and not value.is_integer():
        return f"{int(round(2 * value))}/2"
    return str(int(value)) if isinstance(value, float) else str(value)


def _short_level_label(level):
    qn = level.qn
    if {"N", "J", "F"}.issubset(qn):
        return f"N={_half_label(qn['N'])},J={_half_label(qn['J'])},F={_half_label(qn['F'])}"
    return level.label


def cn_n1_0_demo():
    molecule = load_hfs_molecule(CN_DATA, "CN")
    level_rows, _ = compute_hfs_ncrit(
        molecule,
        [CN_DEMO_TK],
        partner=CN_DEMO_PARTNER,
        require_hfs_collisions=True,
    )
    ncrit_by_upper = {row["upper_id"]: row for row in level_rows}

    # Select the CN N=1 -> 0 hyperfine components and sort them by frequency.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("N") == 1
        and molecule.levels[transition.lower_id].qn.get("N") == 0
    ]
    transitions.sort(key=lambda item: item.frequency_GHz)

    print(f"\nCN N=1-0 HFS components at Tk={CN_DEMO_TK:g} K, partner={CN_DEMO_PARTNER}")
    print("freq_GHz  upper_id  lower_id  upper_level                  lower_level                  ncrit_cm-3")
    for transition in transitions:
        upper = molecule.levels[transition.upper_id]
        lower = molecule.levels[transition.lower_id]
        ncrit = ncrit_by_upper[transition.upper_id]["ncrit_cm-3"]
        print(
            f"{transition.frequency_GHz:8.5f}  "
            f"{transition.upper_id:8d}  "
            f"{transition.lower_id:8d}  "
            f"{_short_level_label(upper):28s} "
            f"{_short_level_label(lower):28s} "
            f"{ncrit:10.3e}"
        )

    _plot_cn_n1_0(transitions, molecule, ncrit_by_upper)
    _plot_cn_n1_0_lte_spectrum(transitions, molecule, ncrit_by_upper)


def _plot_cn_n1_0(transitions, molecule, ncrit_by_upper):
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
        labels.append(
            f"{transition.frequency_GHz:.3f} GHz\n"
            f"u:{_short_level_label(upper)}\n"
            f"l:{_short_level_label(lower)}"
        )
        values.append(ncrit_by_upper[transition.upper_id]["ncrit_cm-3"])

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(values))
    bars = ax.bar(x, values, color="#4C78A8", edgecolor="#243B53", linewidth=0.7)
    ax.set_yscale("log")
    ax.set_ylabel(r"$n_{\rm crit}$ (cm$^{-3}$)")
    ax.set_xlabel("CN N=1-0 hyperfine transition, sorted by frequency")
    ax.set_title(f"CN N=1-0 HFS level-based critical densities at Tk={CN_DEMO_TK:g} K ({CN_DEMO_PARTNER})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=65, ha="right", fontsize=8)
    ax.grid(axis="y", which="both", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.05, f"{value:.1e}", ha="center", va="bottom", fontsize=8, rotation=90)
    fig.tight_layout()
    out_path = OUTPUT_DIR / "cn_n1_0_ncrit_hist.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved CN N=1-0 ncrit histogram: {out_path}")


def _lte_component_weight(transition, molecule):
    upper = molecule.levels[transition.upper_id]
    # Optically thin LTE integrated intensity is proportional to upper population times A_ul.
    return upper.g * np.exp(-upper.energy_K / LTE_TEX) * transition.A_s


def _frequency_fwhm_GHz(frequency_GHz):
    return frequency_GHz * GAUSSIAN_FWHM_KMS / C_KMS


def _split_by_largest_frequency_gap(transitions):
    if len(transitions) < 2:
        return [transitions]
    frequencies = np.array([transition.frequency_GHz for transition in transitions])
    gaps = np.diff(frequencies)
    split_at = int(np.argmax(gaps)) + 1
    if gaps[split_at - 1] > 5.0 * np.median(gaps[gaps > 0]):
        return [transitions[:split_at], transitions[split_at:]]
    return [transitions]


def _plot_cn_n1_0_lte_spectrum(transitions, molecule, ncrit_by_upper):
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
    weights = np.array([_lte_component_weight(transition, molecule) for transition in transitions])
    weights = weights / weights.max()
    weight_by_id = {id(transition): weight for transition, weight in zip(transitions, weights)}
    groups = _split_by_largest_frequency_gap(transitions)
    upper_ids = sorted({transition.upper_id for transition in transitions})
    colors = dict(zip(upper_ids, plt.cm.tab10(np.linspace(0, 1, max(len(upper_ids), 2)))[: len(upper_ids)]))

    group_profiles = []
    norm = 0.0
    for group in groups:
        centers = np.array([transition.frequency_GHz for transition in group])
        max_fwhm = max(_frequency_fwhm_GHz(freq) for freq in centers)
        frequency = np.linspace(centers.min() - 5.0 * max_fwhm, centers.max() + 5.0 * max_fwhm, 2500)
        total = np.zeros_like(frequency)
        components = []
        for transition in group:
            center = transition.frequency_GHz
            sigma = _frequency_fwhm_GHz(center) / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            profile = weight_by_id[id(transition)] * np.exp(-0.5 * ((frequency - center) / sigma) ** 2)
            components.append((transition, center, profile))
            total += profile
        norm = max(norm, total.max())
        group_profiles.append((frequency, total, components))

    fig, axes = plt.subplots(1, len(groups), figsize=(14, 7), sharey=True, constrained_layout=True)
    if len(groups) == 1:
        axes = [axes]
    for ax, (frequency, total, components) in zip(axes, group_profiles):
        for transition, center, profile in components:
            color = colors[transition.upper_id]
            ax.fill_between(frequency, profile / norm, color=color, alpha=0.28)
            ax.plot(frequency, profile / norm, color=color, lw=1.0, alpha=0.85)
            ax.axvline(center, color=color, lw=0.8, alpha=0.65)
        ax.plot(frequency, total / norm, color="black", lw=2.0, label="LTE optically thin sum")
        ax.set_xlabel("Frequency (GHz)")
        ax.grid(alpha=0.25)

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
                label=f"u={_short_level_label(level)}; ncrit={ncrit:.2e} cm^-3",
            )
        )
    axes[-1].legend(handles=handles, fontsize=8, loc="upper right", frameon=False)
    axes[0].set_ylabel("Normalized intensity")
    fig.suptitle(
        f"CN N=1-0 Gaussian LTE optically thin HFS spectrum\n"
        f"Tex={LTE_TEX:g} K, FWHM={GAUSSIAN_FWHM_KMS:g} km s$^{{-1}}$, ncrit at Tk={CN_DEMO_TK:g} K"
    )
    out_path = OUTPUT_DIR / "cn_n1_0_lte_spectrum_ncrit.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved CN N=1-0 LTE HFS spectrum: {out_path}")


if __name__ == "__main__":
    test_branching_and_collision_sums()
    test_detailed_balance_reconstructs_upward_rate()
    cn_n1_0_demo()
    print("toy HFS tests passed")
