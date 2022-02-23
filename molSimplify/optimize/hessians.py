import os
import subprocess
import tempfile
import numpy as np
import ase.io
import ase.units
from molSimplify.optimize.calculators import (_xtb_methods,
                                              _openbabel_methods,
                                              get_calculator)


def compute_guess_hessian(atoms, method):
    if method.lower() in _xtb_methods:
        return xtb_hessian(atoms, method)
    elif method.lower() in _openbabel_methods:
        old_calc = atoms.calc
        atoms.calc = get_calculator(method.lower())
        H = numerical_hessian(atoms)
        atoms.calc = old_calc
        return H


def xtb_hessian(atoms, method):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write .xyz file
        ase.io.write(os.path.join(tmpdir, 'tmp.xyz'), atoms, plain=True)
        try:
            output = subprocess.run(
                ['xtb', '--hess', 'tmp.xyz'],
                cwd=tmpdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        except FileNotFoundError:
            raise ChildProcessError('Could not find subprocess xtb. Ensure xtb'
                                    ' is installed and properly configured.')
        if output.returncode != 0:
            print(output)
            raise ChildProcessError('XTB calculation failed')
        H = read_xtb_hessian(os.path.join(tmpdir, 'hessian'))
    return H


def read_xtb_hessian(file):
    with open(file, 'r') as fin:
        content = fin.read()
    values = np.array([float(f) for f in content.split()[1:]])
    N = int(np.sqrt(values.size))
    return values.reshape(N, N) * ase.units.Hartree / ase.units.Bohr**2


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
