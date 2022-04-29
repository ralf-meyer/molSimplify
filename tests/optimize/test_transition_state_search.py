import pytest
import numpy as np
import ase.build
from molSimplify.optimize.connectivity import get_primitives
from molSimplify.optimize.primitives import (Distance, Angle, Improper,
                                             Dihedral)
from molSimplify.optimize.coordinate_sets import (InternalCoordinates,
                                                  DelocalizedCoordinates)
from molSimplify.optimize.optimizers import RFO
from molSimplify.optimize.hessian_guess import numerical_hessian

from xtb.ase.calculator import XTB


def test_ammonia_transition_state():
    atoms = ase.build.molecule('NH3')
    atoms.calc = XTB(method='GFN2-xTB')

    # The plane spanned by the H atoms is at z=-.27
    atoms.positions[0, :] = [0., 0., -0.1]

    primitives = [Distance(0, 1), Distance(0, 2), Distance(0, 3),
                  Angle(1, 0, 2), Angle(2, 0, 3), Improper(0, 1, 2, 3)]
    H = numerical_hessian(atoms)
    coord_set = InternalCoordinates(primitives)
    opt = RFO(atoms, coordinate_set=coord_set, H0=H, mu=1, maxstep=0.05)
    opt.run(fmax=0.005, steps=100)
    assert opt.converged()
    # Assert that the geometry is close to planar
    assert np.abs(primitives[-1].value(atoms.get_positions())) < 1e-2


def test_ethane_transition_state():
    atoms = ase.build.molecule('C2H6')
    atoms.calc = XTB(method='GFN2-xTB')

    # Rotate one methyl group by 30 degrees about the z axis
    rot_mat = np.array([[np.sqrt(3)/2, -0.5, 0.],
                        [0.5, np.sqrt(3)/2, 0.],
                        [0., 0., 1.]])
    atoms.positions[2:5] = atoms.positions[2:5].dot(rot_mat)

    xyzs = atoms.get_positions()
    bonds = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5), (1, 6), (1, 7)]
    primitives = get_primitives(xyzs, bonds)
    H = numerical_hessian(atoms)
    coord_set = DelocalizedCoordinates(primitives, xyzs)
    opt = RFO(atoms, coordinate_set=coord_set, H0=H, mu=1, maxstep=0.05)
    opt.run(fmax=0.005, steps=100)

    assert opt.converged()
    # Assert that the two methly groups have aligned
    assert abs(Dihedral(3, 0, 1, 6).value(atoms.get_positions())) < 1e-2
    assert abs(Dihedral(2, 0, 1, 7).value(atoms.get_positions())) < 1e-2
    assert abs(Dihedral(4, 0, 1, 5).value(atoms.get_positions())) < 1e-2
