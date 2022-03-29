import os
import pytest
import ase.atoms
import ase.build
import ase.units
import numpy as np
import numdifftools as nd
import geometric.internal
from xtb.ase.calculator import XTB
from molSimplify.optimize.calculators import (_openbabel_methods,
                                              get_calculator)
from molSimplify.optimize.hessians import (filter_hessian,
                                           compute_hessian_guess,
                                           numerical_hessian,
                                           schlegel_hessian)


def num_hessian(atoms, step=None):
    x0 = atoms.get_positions()

    def fun(x):
        atoms.set_positions(x.reshape(-1, 3))
        return -atoms.get_forces().flatten()
    H = nd.Jacobian(fun, step=step)(x0.flatten())
    return 0.5*(H + H.T)


@pytest.mark.parametrize('method', _openbabel_methods)
@pytest.mark.parametrize('system', ['H2', 'H2O', 'C3H8'])
def test_numerical_hessian(method, system):
    atoms = ase.build.molecule(system)
    atoms.calc = get_calculator(method)
    if method == 'mmff94' and system == 'H2':
        # MMFF94 does not have parameters for H2 and is
        # therefore expected to fail.
        with pytest.raises(RuntimeError):
            atoms.get_potential_energy()
    else:
        x0 = atoms.get_positions()
        H = numerical_hessian(atoms, symmetrize=False)
        np.testing.assert_allclose(atoms.get_positions(), x0)
        np.testing.assert_allclose(H, H.T, atol=1e-8)
        H_ref = num_hessian(atoms, step=1e-5)
        # Symmetrize
        H = 0.5*(H + H.T)
        np.testing.assert_allclose(H, H_ref, atol=1e-4)


@pytest.mark.parametrize('system', ['H2', 'LiF'])
def test_xtb_hessian(system):
    """TODO: For some reason I can not get this test to pass for any
    non-dimer molecules. Maybe due to the way the singlepoint calculations are
    restarted in xtb internal hessian calculation. RM 2022/02/22
    """
    atoms = ase.build.molecule(system)
    atoms.rotate(15, (1, 0, 0))
    x0 = atoms.get_positions()
    H = compute_hessian_guess(atoms, 'xtb')
    np.testing.assert_allclose(atoms.get_positions(), x0)
    np.testing.assert_allclose(H, H.T, atol=1e-8)
    atoms.calc = XTB(method='GFN2-xTB', accuracy=0.3)
    H_ref = numerical_hessian(atoms, step=0.005 * ase.units.Bohr)
    eig, _ = np.linalg.eigh(H)
    eig_ref, _ = np.linalg.eigh(H_ref)
    np.testing.assert_allclose(eig, eig_ref, atol=1e-2)
    np.testing.assert_allclose(H, H_ref, atol=1e-2)


@pytest.mark.parametrize('method', ['uff', 'xtb'])
def _test_Fe_CO_6(method):
    """TODO: Fails for both uff (exploding energies) and xtb
    (see test_xtb_hessian)
    """
    r_FeC = 2.3
    r_FeO = 2.3 + 1.1
    atoms = ase.atoms.Atoms(['Fe']+['C', 'O']*6,
                            positions=[[0., 0., 0.],
                                       [r_FeC, 0., 0.], [r_FeO, 0., 0.],
                                       [0., r_FeC, 0.], [0., r_FeO, 0.],
                                       [-r_FeC, 0., 0.], [-r_FeO, 0., 0.],
                                       [0., -r_FeC, 0.], [0., -r_FeO, 0.],
                                       [0., 0., r_FeC], [0., 0., r_FeO],
                                       [0., 0., -r_FeC], [0., 0., -r_FeO]],
                            charges=[2]+[0, 0]*6)
    x0 = atoms.get_positions()
    atoms.calc = get_calculator(method)
    H = compute_hessian_guess(atoms, method)
    np.testing.assert_allclose(atoms.get_positions(), x0)
    np.testing.assert_allclose(H, H.T, atol=1e-8)
    H_ref = num_hessian(atoms, step=1e-5)
    np.testing.assert_allclose(H, H_ref, atol=1e-4)


@pytest.mark.parametrize('system', ['H2', 'F2', 'H2O', 'CO2'])
def test_schlegel_vs_geometric(tmpdir, system):
    """ Only tests simple systems since geomeTRIC does not use the actual
    Schlegel rules for dihedrals and impropers."""
    atoms = ase.build.molecule(system)
    tmp_file = os.path.join(tmpdir, 'tmp.xyz')
    ase.io.write(tmp_file, atoms, plain=True)
    mol = geometric.molecule.Molecule(tmp_file)
    coords_ref = geometric.internal.PrimitiveInternalCoordinates(
        mol, connect=True)
    xyzs = atoms.get_positions()
    # The following only works in the newest version of geometric:
    # H_ref = coords_ref.calcHessCart(
    #     xyzs/ase.units.Bohr, np.zeros(len(coords_ref.Internals)),
    #     coords_ref.guess_hessian(xyzs/ase.units.Bohr))
    # for now calcHessCart is replaced with a simplified function
    # that neglects the gradient in internal coordinates:
    Bmat = coords_ref.wilsonB(xyzs/ase.units.Bohr)
    H_ref = np.einsum('ai,ab,bj->ij', Bmat,
                      coords_ref.guess_hessian(xyzs/ase.units.Bohr),
                      Bmat) * ase.units.Hartree/ase.units.Bohr**2
    H = schlegel_hessian(atoms, threshold=1.2)
    np.testing.assert_allclose(H, H_ref, atol=1e-10)


def test_filter_hessian():
    H = np.diag([-1., 0., 1., 2., 3., 4.])
    H = filter_hessian(H, thresh=1.1e-5)
    np.testing.assert_allclose(H, np.diag([1.1e-5, 1.1e-5, 1., 2., 3., 4.]))

    # Build random matrix with eigenvalue above 0.1
    A = np.array([[0.7432, 0.4965, 0.2700, 0.0742, -0.0800, -0.1814],
                  [0.4965, 1.0133, 0.5708, 0.1900, -0.1071, -0.2977],
                  [0.2700, 0.5708, 0.9332, 0.3893, -0.0277, -0.2830],
                  [0.0742, 0.1900, 0.3893, 0.7155,  0.2134, -0.0696],
                  [-0.0800, -0.1071, -0.0277, 0.2134, 0.6736, 0.4129],
                  [-0.1814, -0.2977, -0.2830, -0.0696, 0.4129, 1.2388]])
    A_ref = A.copy()
    A = filter_hessian(A)
    # Test that it remained unaltered
    np.testing.assert_allclose(A, A_ref)
    # Increase the threshold to 1.0
    A = filter_hessian(A, thresh=1.0)
    vals, _ = np.linalg.eigh(A)
    # Test that the smaller eigenvalues have been filtered correctly
    np.testing.assert_allclose(vals, [1., 1., 1., 1., 1.27649882, 2.20586986])
