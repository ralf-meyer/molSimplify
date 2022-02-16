import ase.io
import ase.constraints
import pytest
import numpy as np
from molSimplify.optimize.optimize import run_preoptimization
from molSimplify.optimize.calculators import _available_calculators
from pkg_resources import resource_filename, Requirement


@pytest.mark.parametrize('method', _available_calculators)
def test_acac(method):
    in_file = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/optimize/inputs/acac/fe_oct_2_acac_3_s_5_conf_1.xyz')
    ref_file = resource_filename(Requirement.parse('molSimplify'),
                                 f'tests/optimize/ref/acac/{method}_opt.xyz')
    atoms = ase.io.read(in_file)
    atoms.set_constraint(ase.constraints.FixAtoms(
        indices=[0, 1, 6, 15, 20, 29, 34]))
    # Assign spin and charge to first atom
    q = np.zeros_like(atoms.get_initial_charges())
    q[0] = -1.0
    atoms.set_initial_charges(q)
    s = np.zeros_like(atoms.get_initial_magnetic_moments())
    # Number of unpaired electrons is multiplicity - 1
    s[0] = 4
    atoms.set_initial_magnetic_moments(s)
    if method == 'mmff94':
        # MMFF94 does not have parameters for Fe and is
        # therefore expected to fail.
        with pytest.raises(RuntimeError):
            run_preoptimization(atoms, method)
    else:
        run_preoptimization(atoms, method)
    atoms_ref = ase.io.read(ref_file)
    ase.build.minimize_rotation_and_translation(atoms, atoms_ref)
    rmsd = np.sqrt(np.mean((atoms.positions - atoms_ref.positions)**2))
    assert rmsd < 0.1
