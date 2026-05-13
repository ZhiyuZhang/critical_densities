"""HCN HFS critical-density demo using hcn@hfs.dat."""

import os
from pathlib import Path

import numpy as np

from critical_density_hfs import compute_hfs_ncrit, load_hfs_molecule


BASE_DIR = Path(__file__).resolve().parent
DATA_FILE = BASE_DIR / "hcn@hfs.dat"
OUTPUT_DIR = BASE_DIR / "outputs"
DEMO_TK = 30.0
DEMO_PARTNER = "p-H2"


def _value_label(value):
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _level_label(level):
    qn = level.qn
    if {"J", "F"}.issubset(qn):
        return f"J={_value_label(qn['J'])},F={_value_label(qn['F'])}"
    return level.label


def _select_hcn_j1_0(molecule):
    # HCN J=1-0 HFS components are all radiative branches with upper J=1 and lower J=0.
    transitions = [
        transition
        for transition in molecule.radiative
        if molecule.levels[transition.upper_id].qn.get("J") == 1
        and molecule.levels[transition.lower_id].qn.get("J") == 0
    ]
    return sorted(transitions, key=lambda item: item.frequency_GHz)


def _compute_rows(molecule):
    level_rows, _ = compute_hfs_ncrit(
        molecule,
        [DEMO_TK],
        partner=DEMO_PARTNER,
        require_hfs_collisions=True,
    )
    return {row["upper_id"]: row for row in level_rows}


def _print_table(transitions, molecule, ncrit_by_upper):
    print(f"\nHCN J=1-0 HFS components at Tk={DEMO_TK:g} K, partner={DEMO_PARTNER}")
    print("freq_GHz  upper_id  lower_id  upper_level        lower_level        ncrit_cm-3")
    for transition in transitions:
        upper = molecule.levels[transition.upper_id]
        lower = molecule.levels[transition.lower_id]
        ncrit = ncrit_by_upper[transition.upper_id]["ncrit_cm-3"]
        print(
            f"{transition.frequency_GHz:8.5f}  "
            f"{transition.upper_id:8d}  "
            f"{transition.lower_id:8d}  "
            f"{_level_label(upper):18s} "
            f"{_level_label(lower):18s} "
            f"{ncrit:10.3e}"
        )


def _plot(transitions, molecule, ncrit_by_upper):
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

    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(values))
    bars = ax.bar(x, values, color="#59A14F", edgecolor="#244A22", linewidth=0.7)
    ax.set_ylabel(r"$n_{\rm crit}$ (cm$^{-3}$)")
    ax.set_xlabel("HCN J=1-0 hyperfine transition, sorted by frequency")
    ax.set_title(f"HCN J=1-0 HFS level-based critical densities at Tk={DEMO_TK:g} K ({DEMO_PARTNER})")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.01, f"{value:.1e}", ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    out_path = OUTPUT_DIR / "hcn_j1_0_ncrit_hist.png"
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"\nSaved HCN J=1-0 ncrit histogram: {out_path}")


def main():
    molecule = load_hfs_molecule(DATA_FILE, "HCN")
    transitions = _select_hcn_j1_0(molecule)
    ncrit_by_upper = _compute_rows(molecule)
    _print_table(transitions, molecule, ncrit_by_upper)
    _plot(transitions, molecule, ncrit_by_upper)


if __name__ == "__main__":
    main()
