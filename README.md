# critical_densities

Calculate optically thin, no-background critical densities following Shirley (2015, PASP, 127, 299):

https://arxiv.org/abs/1501.01629

The implemented quantity is the Shirley equation (4) multi-level critical density:

```text
n_crit = A_jk / sum_i gamma_ji
```

where the denominator includes all collisional depopulation routes out of the upper level `j`. Downward rates are read directly from the LAMDA data files. Upward rates are computed from downward rates with detailed balance, Shirley equation (5):

```text
gamma_ji = gamma_ij * (g_i / g_j) * exp(-(E_i - E_j) / T_kin)
```

This code does not calculate `n_eff`, optical-depth corrections, radiative trapping, CMB/background radiation corrections, or dust/IR pumping corrections.

## Data

The calculations use local LAMDA-style data files in this directory:

```text
12co.dat
13co.dat
c18o.dat
hcn.dat
hco+.dat
hnc.dat
cs.dat
cn.dat
```

The parser now reads the LAMDA sections by header text instead of hard-coded row numbers, so small formatting differences between files are handled more safely.

## Usage

From the repository root:

```python
import sys
sys.path.insert(0, "critical_densities")

from critical_density import ncrit

print(ncrit("HCN", 0, 10, 3, False))
print(ncrit("HCO+", 2, 50, 3, True))
```

The legacy call signature is still supported:

```python
ncrit(molecule, J_low, Tkin, op_ratio, verbose)
```

For simple linear rotors, `J_low=0` means `J=1-0`, `J_low=1` means `J=2-1`, and so on.
Verbose output reports both the quantum-number transition and the underlying LAMDA level indices; for example, `2->1 (LAMDA levels 3->2)` is the `J=2-1` line.

For spectra that are not a simple one-level-per-J ladder, such as CN, select the LAMDA radiative transition explicitly:

```python
ncrit("CN", T=10, transition_index=2)
ncrit("CN", T=10, upper=3, lower=1)
```

## Fixes in This Version

- Interpolation is now applied to both downward and upward collision terms.
- Non-grid temperatures such as 15 K and 25 K no longer drop the downward collision sum.
- Ortho/para H2 mixing is now done on the collision-rate denominator before dividing by `A_jk`.
- Data files are found relative to `critical_density.py`, so imports work from the repository root.
- The old debug `print(S)` output has been removed.
- The upward-collision sum now uses the actual levels in the file instead of stopping at a hard-coded level index.
- HNC, CS, and CN are supported through the general LAMDA parser.
- Electron collision blocks are ignored by default because Shirley Table 1 reports collisions with H2 or H2-scaled neutral partners.

## Notes and Limits

Temperatures must lie within the tabulated range of the selected neutral collision partner. The code interpolates but does not silently extrapolate.

Some Shirley table entries use extrapolated rates when the LAMDA grid does not include a requested temperature. This implementation raises a clear error instead of extrapolating unless the data file itself covers that temperature.

The `op_ratio` argument only matters when separate ortho-H2 and para-H2 collision partners are present, as in the CO isotopologue files. For files with a single H2, He-scaled-H2, or He-like neutral partner, that partner is used with weight 1.

## HFS-Resolved Critical Density

For hyperfine-resolved files such as `cn-hfs.dat`, use the separate script:

```bash
python critical_density_hfs.py \
    --data cn-hfs.dat \
    --species CN \
    --partner p-H2 \
    --temperatures 5,10,15,20,25,30,35,40,45,50 \
    --require-hfs-collisions \
    --out outputs \
    --verbose
```

The HFS script computes a level-based critical density:

```text
ncrit(u; Tk) = A_tot(u) / gamma_tot(u; Tk)
A_tot(u) = sum_l A_ul
gamma_tot(u; Tk) = sum_i gamma_ui(Tk)
```

Here `u` is an HFS level such as CN `(N,J,F)`. Multiple satellite components with the same upper level share the same `ncrit`; their individual `A_ul` values only define branching ratios.

The script writes:

```text
level_ncrit.csv
transition_branches.csv
```

`level_ncrit.csv` has one row per upper HFS level and temperature. It includes `upper_group`, for example `N=1`, `N=2`, and `N=3` for CN, so HFS sublevels remain easy to scan by parent rotational level. `transition_branches.csv` has one row per radiative component and includes `branching_ratio = A_ul / A_tot`.

The script does not silently mix para-H2 and ortho-H2. Select one partner with `--partner p-H2` or `--partner o-H2`, or explicitly request a mixture:

```bash
python critical_density_hfs.py \
    --data cn-hfs.dat \
    --species CN \
    --h2-opr 3 \
    --temperatures 5,10,15,20,25,30,35,40,45,50 \
    --require-hfs-collisions \
    --out outputs
```

Temperature interpolation uses log-log interpolation for strictly positive collision rates. Temperatures outside the tabulated range raise an error unless `--allow-extrapolation` is supplied, in which case the CSV warnings field marks the extrapolation.

### HFS Demo Scripts

Small plotting demos are available for CN, HCN, and NH3:

```bash
python cn_hfs_toy.py
python hcn_hfs_toy.py
python nh3_hfs_toy.py
```

The scripts print the selected HFS components and save bar charts under `outputs/`.

- `cn_hfs_toy.py`: CN `N=1-0`, using `cn-hfs.dat`.
- `hcn_hfs_toy.py`: HCN `J=1-0`, using `hcn@hfs.dat`; default `Tk=30 K` because that collision table currently spans 5-30 K.
- `nh3_hfs_toy.py`: p-NH3 `(1,1)` inversion and o-NH3 `1_0-0_0`, using the Loreau et al. p/o-NH3 files; default `Tk=50 K`.
