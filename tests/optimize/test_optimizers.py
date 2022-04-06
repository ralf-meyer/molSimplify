import pytest
import numpy as np
import ase.atoms
import ase.build
import ase.calculators.emt
import ase.optimize
from molSimplify.optimize.calculators import OpenbabelFF
from molSimplify.optimize.connectivity import find_connectivity, get_primitives
from molSimplify.optimize.primitives import Distance
from molSimplify.optimize.coordinate_sets import (DelocalizedCoordinates,
                                                  InternalCoordinates)
from molSimplify.optimize.optimizers import NonCartesianBFGS, NonCartesianLBFGS
from molSimplify.Scripts.rmsd import kabsch_rmsd


@pytest.mark.parametrize('opt', [NonCartesianBFGS, NonCartesianLBFGS])
def test_optimizers_on_H2(opt):
    """Separate test case because H2 is not supported by MMFF94"""
    atoms = ase.atoms.Atoms(['H', 'H'], positions=[[0., 0., 0.],
                                                   [.5, .5, .5]])
    atoms.calc = ase.calculators.emt.EMT()

    # Reference calculation in Cartesian coordinates
    atoms_ref = atoms.copy()
    atoms_ref.calc = ase.calculators.emt.EMT()
    opt_ref = ase.optimize.LBFGS(atoms_ref)
    opt_ref.run(fmax=0.01)
    xyzs_ref = atoms_ref.get_positions()
    r_ref = np.linalg.norm(xyzs_ref[0] - xyzs_ref[1])

    coord_set = InternalCoordinates([Distance(0, 1)])
    opt = NonCartesianLBFGS(atoms, coord_set)
    opt.run(fmax=0.01)
    # Test that the final bond length correct
    xyzs = atoms.get_positions()
    r = np.linalg.norm(xyzs[0] - xyzs[1])
    assert abs(r - r_ref) < 1e-3


@pytest.mark.parametrize('opt, mol, coord_set',
                         [(opt, mol, coord_set)
                          for opt in ['BFGS', 'LBFGS']
                          for mol in ['H2O', 'NH3', 'CH4', 'C2H4', 'C2H6',
                                      'C6H6', 'butadiene', 'bicyclobutane']
                          for coord_set in ['internal', 'dlc']])
def test_optimizers_on_organic_molecules(opt, mol, coord_set):
    atoms = ase.build.molecule(mol)
    # slightly distort
    xyzs = atoms.get_positions()
    rng = np.random.default_rng(1234)
    xyzs += rng.normal(scale=0.05, size=xyzs.shape)
    atoms.set_positions(xyzs)
    atoms.calc = OpenbabelFF(ff='MMFF94')

    # Reference calculation in Cartesian coordinates
    atoms_ref = atoms.copy()
    atoms_ref.calc = OpenbabelFF(ff='MMFF94')

    bonds = find_connectivity(atoms)
    primitives = get_primitives(xyzs, bonds)
    if coord_set == 'internal':
        coord_set = InternalCoordinates(primitives)
    elif coord_set == 'dlc':
        coord_set = DelocalizedCoordinates(primitives, xyzs)

    if opt == 'BFGS':
        opt_ref = ase.optimize.BFGS(atoms_ref)
        opt = NonCartesianBFGS(atoms, coord_set)
    elif opt == 'LBFGS':
        opt_ref = ase.optimize.LBFGS(atoms_ref)
        opt = NonCartesianLBFGS(atoms, coord_set)

    opt_ref.run(fmax=0.001, steps=300)
    opt.run(fmax=0.001, steps=300)

    assert kabsch_rmsd(atoms.get_positions(),
                       atoms_ref.get_positions(),
                       translate=True) < 1e-2
