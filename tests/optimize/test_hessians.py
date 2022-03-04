import pytest
import ase.atoms
import ase.build
import numpy as np
import numdifftools as nd
from molSimplify.optimize.calculators import (_available_calculators,
                                              get_calculator)
from molSimplify.optimize.hessians import filter_hessian, numerical_hessian


def ref_hessian(atoms, step=None):
    x0 = atoms.get_positions()

    def fun(x):
        atoms.set_positions(x.reshape(-1, 3))
        return -atoms.get_forces().flatten()
    H = nd.Jacobian(fun, step=step)(x0.flatten())
    return 0.5*(H + H.T)


@pytest.mark.parametrize('method', _available_calculators)
def test_numerical_hessian_H2(method):
    h2 = ase.build.molecule('H2')
    h2.calc = get_calculator(method)
    if method == 'mmff94':
        # MMFF94 does not have parameters for H2 and is
        # therefore expected to fail.
        with pytest.raises(RuntimeError):
            h2.get_potential_energy()
    else:
        x0 = h2.get_positions()
        H = numerical_hessian(h2, symmetrize=False)
        np.testing.assert_allclose(h2.get_positions(), x0)
        np.testing.assert_allclose(H, H.T, atol=1e-8)
        H_ref = ref_hessian(h2)
        # Symmetrize
        H = 0.5*(H + H.T)
        np.testing.assert_allclose(H, H_ref, atol=1e-4)


@pytest.mark.parametrize('method', _available_calculators)
def test_numerical_hessian_water(method):
    atoms = ase.build.molecule('H2O')
    atoms.calc = get_calculator(method)
    x0 = atoms.get_positions()
    H = numerical_hessian(atoms, symmetrize=False)
    np.testing.assert_allclose(atoms.get_positions(), x0)
    np.testing.assert_allclose(H, H.T, atol=1e-8)
    H_ref = ref_hessian(atoms)
    # Symmetrize
    H = 0.5*(H + H.T)
    np.testing.assert_allclose(H, H_ref, atol=1e-4)


@pytest.mark.parametrize('method', _available_calculators)
def test_numerical_hessian_benzene(method):
    atoms = ase.build.molecule('C6H6')
    atoms.calc = get_calculator(method)
    x0 = atoms.get_positions()
    H = numerical_hessian(atoms, symmetrize=False)
    np.testing.assert_allclose(atoms.get_positions(), x0)
    np.testing.assert_allclose(H, H.T, atol=1e-8)
    H_ref = ref_hessian(atoms, step=1e-5)
    # Symmetrize
    H = 0.5*(H + H.T)
    np.testing.assert_allclose(H, H_ref, atol=1e-4)


@pytest.mark.parametrize('method', _available_calculators)
def test_numerical_hessian_Fe_CO_6(method):
    atoms = ase.atoms.Atoms(['Fe']+['C', 'O']*6,
                            positions=[[0., 0., 0.],
                                       [2.3, 0., 0.], [3.4, 0., 0.],
                                       [0., 2.3, 0.], [0., 3.4, 0.],
                                       [-2.3, 0., 0.], [-3.4, 0., 0.],
                                       [0., -2.3, 0.], [0., -3.4, 0.],
                                       [0., 0., 2.3], [0., 0., 3.4],
                                       [0., 0., -2.3], [0., 0., -3.4]],
                            charges=[2]+[0, 0]*6)
    atoms.calc = get_calculator(method)
    if method == 'mmff94':
        # MMFF94 does not have parameters for Fe and is
        # therefore expected to fail.
        with pytest.raises(RuntimeError):
            atoms.get_potential_energy()
    else:
        x0 = atoms.get_positions()
        H = numerical_hessian(atoms, symmetrize=False)
        np.testing.assert_allclose(atoms.get_positions(), x0)
        np.testing.assert_allclose(H, H.T, atol=1e-8)
        H_ref = ref_hessian(atoms, step=1e-5)
        # Symmetrize
        H = 0.5*(H + H.T)
        np.testing.assert_allclose(H, H_ref, atol=1e-4)


def test_filter_hessian():
    H = np.diag([-1., 0., 1., 2., 3., 4.])
    H = filter_hessian(H)
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
