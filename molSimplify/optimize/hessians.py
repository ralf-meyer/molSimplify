import os
import subprocess
import tempfile
import numpy as np
import ase.io
import ase.units
from molSimplify.optimize.calculators import (_xtb_methods,
                                              _openbabel_methods,
                                              get_calculator)
from molSimplify.optimize.connectivity import (find_connectivity,
                                               find_primitives)


def compute_guess_hessian(atoms, method):
    if method.lower() in _xtb_methods:
        return xtb_hessian(atoms, method)
    elif method.lower() in _openbabel_methods:
        old_calc = atoms.calc
        atoms.calc = get_calculator(method.lower())
        H = numerical_hessian(atoms)
        atoms.calc = old_calc
        return H


def schlegel_hessian(atoms):
    """
    Schlegel, Theoret. Chim. Acta 66, 333-340 (1984).
    https://doi.org/10.1007/BF00554788

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        Arrangement of atoms.
    Returns
    -------
    H : np.ndarray
        Guess Hessian in cartesian coordinates and ase units (eV, Ang)
    """
    atomic_numbers = atoms.get_atomic_numbers()
    xyzs = atoms.get_positions()
    # Calculate the covalent bond distances as they are needed later
    cov = np.array([ase.data.covalent_radii[num] for num in atomic_numbers])
    r_cov = cov[:, np.newaxis] + cov[np.newaxis, :]
    # "Atoms are considered bonded if their internuclear distance is less than
    # 1.35 times the sum of the covalentradii"
    bonds = find_connectivity(atoms, threshold=1.35**2, connect_fragments=True)
    bends, linear_bends, torsions, planars = find_primitives(xyzs, bonds)

    def get_B(num1, num2):
        """Returns the B parameter given two atomic numbers"""
        # Sort for simplicity
        num1, num2 = min(num1, num2), max(num1, num2)
        if num1 <= 2:  # First period
            if num2 <= 2:  # first-first
                return -0.244
            elif num2 <= 10:  # first-second
                return 0.352
            else:  # first-third+
                return 0.660
        elif num1 <= 10:  # Second period
            if num2 <= 10:  # second-second
                return 1.085
            else:  # second-third+
                return 1.522
        else:  # third+-third+
            return 2.068

    F_str = []
    for b in bonds:
        r = np.linalg.norm(xyzs[b.i] - xyzs[b.k])
        if r < 1.35 * r_cov[b.i, b.k]:
            B = get_B(atomic_numbers[b.i], atomic_numbers[b.j])
            F_str.append(1.734/(r - B)**3)
        else:
            # Not covalently bonded atoms (from fragment connection algorithm).
            # Following geomeTRIC those are assigned a fixed value:
            F_str.append(0.1 * ase.units.Hartree)
    F_bend = []
    for a in bends + linear_bends:
        if atomic_numbers[a.i] == 1 or atomic_numbers[a.k] == 1:
            # "either or both terminal atoms hydrogen"
            F_bend.append(0.160)
        else:
            # "all three heavy atom bends"
            F_bend.append(0.250)
    F_tors = []
    for t in torsions:
        r = np.linalg.norm(xyzs[t.j] - xyzs[t.k])
        F_tors = 0.0023 - 0.07*(r - r_cov[t.j, t.k])
    F_oop = []
    for p in planars:
        r1 = xyzs[p.j] - xyzs[p.i]
        r2 = xyzs[p.k] - xyzs[p.i]
        r3 = xyzs[p.l] - xyzs[p.i]
        # Additional np.abs() since we do not know the orientation of r1
        # with respect to r2 x r3.
        d = 1 - np.abs(np.dot(r1, np.cross(r2, r3)))/(
            np.linalg.norm(r1)*np.linalg.norm(r2)*np.linalg.norm(r3))
        F_oop.append(0.045 * d**4)

    H = np.diag(F_str + F_bend + F_tors + F_oop)
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
