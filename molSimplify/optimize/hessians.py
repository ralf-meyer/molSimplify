import numpy as np
import numdifftools as nd
from molSimplify.optimize.calculators import (_available_calculators,
                                              get_calculator)


def compute_guess_hessian(atoms, guess_hessian):
    if guess_hessian.lower() in _available_calculators:
        atoms.calc = get_calculator(guess_hessian.lower())
        return numerical_hessian(atoms)


def numerical_hessian(atoms, step=1e-5, symmetrize=True):
    N = len(atoms)
    x0 = atoms.get_positions()
    H = np.zeros((3*N, 3*N))

    for i in range(N):
        for c in range(3):
            x = x0.copy()
            x[i, c] += step
            atoms.set_positions(x)
            g_plus = -atoms.get_forces().flatten()

            x = x0.copy()
            x[i, c] -= step
            atoms.set_positions(x)
            g_minus = -atoms.get_forces().flatten()
            H[3*i + c, :] = (g_plus - g_minus)/(2*step)
    atoms.set_positions(x0)
    if symmetrize:
        return 0.5*(H + H.T)
    return H


def filter_hessian(H, thresh=1.1e-5):
    """GeomeTRIC resets calculations if Hessian eigenvalues below
    a threshold of 1e-5 are encountered. This method is used to
    construct a new Hessian matrix where all eigenvalues smaller
    than the threshold are set exactly to the threshold value
    which by default is slightly above geomeTRICs cutoff.

    Parameters
    ----------
    H : np.array
        input Hessian
    thresh : float
        filter threshold. Default 1.1e-5

    Returns
    -------
    H : np.array
        filtered Hessian
    """
    vals, vecs = np.linalg.eigh(H)
    vals[vals < thresh] = thresh
    H = np.einsum('ji,i,ki->jk', vecs, vals, vecs)
    return H
