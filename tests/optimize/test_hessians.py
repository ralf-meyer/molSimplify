import pytest
import ase.atoms
import ase.build
import numpy as np
import numdifftools as nd
from molSimplify.optimize.calculators import (_available_calculators,
                                              get_calculator)
from molSimplify.optimize.hessians import numerical_hessian


def ref_hessian(atoms, step=None):
    x0 = atoms.get_positions()

    def fun(x):
        atoms.set_positions(x.reshape(-1, 3))
        return atoms.get_potential_energy()
    atoms.set_positions(x0)
    return nd.Hessian(fun, step=step)(x0.flatten())


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
    x0 = atoms.get_positions()
    H = numerical_hessian(atoms, symmetrize=False)
    np.testing.assert_allclose(atoms.get_positions(), x0)
    np.testing.assert_allclose(H, H.T, atol=1e-8)
    H_ref = ref_hessian(atoms, step=1e-5)
    # Symmetrize
    H = 0.5*(H + H.T)
    np.testing.assert_allclose(H, H_ref, atol=1e-4)
