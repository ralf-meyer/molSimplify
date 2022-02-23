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
from molSimplify.optimize.coordinates import (Distance, Angle,
                                              LinearAngle, Dihedral)


def compute_guess_hessian(atoms, method):
    if method.lower() in _xtb_methods:
        return xtb_hessian(atoms, method)
    elif method.lower() in _openbabel_methods:
        old_calc = atoms.calc
        atoms.calc = get_calculator(method.lower())
        H = numerical_hessian(atoms)
        atoms.calc = old_calc
        return H
    elif method.lower() == 'schlegel':
        return schlegel_hessian(atoms)


def schlegel_hessian(atoms, threshold=1.35):
    """
    Schlegel, Theoret. Chim. Acta 66, 333-340 (1984).
    https://doi.org/10.1007/BF00554788

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        Arrangement of atoms.
    threshold : float
        Atoms closer than this threshold times the sum of their covalent radii
        are considered bound. Default value suggested by Schlegel is 1.35.

    Returns
    -------
    H : numpy.ndarray
        Guess Hessian in cartesian coordinates and ase units (eV/Ang**2)
    """
    N = len(atoms)
    atomic_numbers = atoms.get_atomic_numbers()
    xyzs = atoms.get_positions()
    # Calculate the covalent bond distances as they are needed later
    cov = np.array([ase.data.covalent_radii[num] for num in atomic_numbers])
    r_cov = cov[:, np.newaxis] + cov[np.newaxis, :]
    # "Atoms are considered bonded if their internuclear distance is less than
    # 1.35 times the sum of the covalentradii"
    bonds = find_connectivity(atoms, threshold=threshold**2,
                              connect_fragments=True)
    bends, linear_bends, torsions, planars = find_primitives(xyzs, bonds)
    # Initialize Hessian in Cartesian coordinates
    H = np.zeros((3*N, 3*N))

    def get_B_str(num1, num2):
        """Returns the B parameter for stretch coordinates given two atomic
        numbers
        """
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

    for b in bonds:
        prim = Distance(*b)
        r = np.linalg.norm(xyzs[b[0]] - xyzs[b[1]])
        if r < threshold * r_cov[b[0], b[1]]:
            B = get_B_str(atomic_numbers[b[0]], atomic_numbers[b[1]])
            h_str = (1.734 * ase.units.Hartree * ase.units.Bohr /
                     (r - B*ase.units.Bohr)**3)
        else:
            # Not covalently bonded atoms (from fragment connection algorithm).
            # Following geomeTRIC those are assigned a fixed value:
            h_str = 0.1 * ase.units.Hartree / ase.units.Bohr**2
        Bi = prim.derivative(xyzs)
        H += np.outer(Bi, h_str * Bi)

    def get_A_bend(num1, num2):
        if atomic_numbers[num1] == 1 or atomic_numbers[num2] == 1:
            # "either or both terminal atoms hydrogen"
            return 0.160 * ase.units.Hartree
        # "all three heavy atom bends"
        return 0.250 * ase.units.Hartree

    for a in bends:
        prim = Angle(*a)
        h_bend = get_A_bend(prim.i, prim.k)
        Bi = prim.derivative(xyzs)
        H += np.outer(Bi, h_bend * Bi)
    for a in linear_bends:
        # Primitives for both projection axes
        prim_u = LinearAngle(*a, axis=0)
        prim_w = LinearAngle(*a, axis=1)
        h_bend = get_A_bend(prim_u.i, prim_u.k)
        Bi = prim_u.derivative(xyzs)
        H += np.outer(Bi, h_bend * Bi)
        Bi = prim_w.derivative(xyzs)
        H += np.outer(Bi, h_bend * Bi)

    for t in torsions:
        prim = Dihedral(*t)
        r = np.linalg.norm(xyzs[t[1]] - xyzs[t[2]])
        h_tors = ase.units.Hartree * (
            0.0023 - 0.07 / ase.units.Bohr * (r - r_cov[t[1], t[2]]))
        Bi = prim.derivative(xyzs)
        H += np.outer(Bi, h_tors * Bi)

    for p in planars:
        prim = Dihedral(*p)
        r1 = xyzs[p[1]] - xyzs[p[0]]
        r2 = xyzs[p[2]] - xyzs[p[0]]
        r3 = xyzs[p[3]] - xyzs[p[0]]
        # Additional np.abs() since we do not know the orientation of r1
        # with respect to r2 x r3.
        d = 1 - np.abs(np.dot(r1, np.cross(r2, r3)))/(
            np.linalg.norm(r1)*np.linalg.norm(r2)*np.linalg.norm(r3))
        h_oop = 0.045 * ase.units.Hartree * d**4
        Bi = prim.derivative(xyzs)
        H += np.outer(Bi, h_oop * Bi)
    return H


def fischer_almloef_hessian(atoms):
    """
    Fischer and Almloef, J. Phys. Chem. 1992, 96, 24, 9768-9774
    https://doi.org/10.1021/j100203a036

    Parameters
    ----------
    atoms : ase.atoms.Atoms
        Arrangement of atoms.

    Returns
    -------
    H : numpy.ndarray
        Guess Hessian in cartesian coordinates and ase units (eV/Ang**2)
    """
    N = len(atoms)
    H = np.zeros((3*N, 3*N))
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
