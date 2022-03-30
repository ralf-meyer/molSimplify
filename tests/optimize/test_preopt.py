import ase.io
import ase.constraints
import pytest
import numpy as np
from molSimplify.optimize.main import run_preoptimization
from molSimplify.optimize.calculators import (_available_methods,
                                              _openbabel_methods)
from pkg_resources import resource_filename, Requirement


@pytest.mark.parametrize('method', _available_methods)
def test_acac(method):
    if method in _openbabel_methods:
        # check if openbabel version > 3.0. This is necessary as
        # OBForceField.GetGradient is not public for prior versions.
        pytest.importorskip('openbabel', minversion='3.0')
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
    run_preoptimization(atoms, method)
    atoms_ref = ase.io.read(ref_file)
    ase.build.minimize_rotation_and_translation(atoms, atoms_ref)
    rmsd = np.sqrt(np.mean((atoms.positions - atoms_ref.positions)**2))
    assert rmsd < 0.1
