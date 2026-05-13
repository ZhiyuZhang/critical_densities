"""Small toy checks for the HFS critical-density implementation."""

import math
from pathlib import Path

import numpy as np

from critical_density_hfs import (
    CollisionPartner,
    HFSMolecule,
    HFSLevel,
    RadiativeTransition,
    _collision_sum_for_level,
    compute_hfs_ncrit,
)


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


if __name__ == "__main__":
    test_branching_and_collision_sums()
    test_detailed_balance_reconstructs_upward_rate()
    print("toy HFS tests passed")
