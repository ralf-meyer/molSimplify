import numpy as np
import ase
from abc import abstractmethod


class Primitive():

    @abstractmethod
    def value(self, xyzs):
        ...

    @abstractmethod
    def derivative(self, xyzs):
        ...


class Distance(Primitive):

    def __init__(self, i, j):
        self.i, self.j = i, j

    def value(self, xyzs):
        rij = xyzs[self.i, :] - xyzs[self.j, :]
        return np.linalg.norm(rij)

    def derivative(self, xyzs):
        rij = xyzs[self.i, :] - xyzs[self.j, :]
        r = np.linalg.norm(rij)
        dr = np.zeros(xyzs.size)
        dr[3*self.i:3*(self.i+1)] = rij/r
        dr[3*self.j:3*(self.j+1)] = -rij/r
        return dr


class InverseDistance(Distance):

    def value(self, xyzs):
        return 1./Distance.q(self, xyzs)

    def derivative(self, xyzs):
        rij = xyzs[self.i, :] - xyzs[self.j, :]
        q = 1.0/np.linalg.norm(rij)
        dq = np.zeros(len(xyzs)*3)
        dq[3*self.i:3*(self.i+1)] = -rij * q**3
        dq[3*self.j:3*(self.j+1)] = rij * q**3
        return dq


class Angle(Primitive):

    def __init__(self, i, j, k):
        """Atom j is center of angle
        Implementation follows:
        https://www.cs.utexas.edu/users/evouga/uploads/4/5/6/8/45689883/turning.pdf
        """
        self.i, self.j, self.k = i, j, k

    def value(self, xyzs):
        rji = xyzs[self.i, :] - xyzs[self.j, :]
        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        norm_rji = np.linalg.norm(rji)
        norm_rjk = np.linalg.norm(rjk)
        cross = np.cross(rji, rjk)
        z = cross/np.linalg.norm(cross)
        t = 2 * np.arctan2(
            np.dot(cross, z),
            norm_rji*norm_rjk + np.dot(rji, rjk))
        return t

    def derivative(self, xyzs):
        rji = xyzs[self.i, :] - xyzs[self.j, :]
        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        norm_rji = np.linalg.norm(rji)
        norm_rjk = np.linalg.norm(rjk)
        cross = np.cross(rji, rjk)
        z = cross/np.linalg.norm(cross)
        dt = np.zeros(xyzs.size)
        dt[3*self.i:3*(self.i+1)] = np.cross(rji, z)/norm_rji**2
        dt[3*self.k:3*(self.k+1)] = - np.cross(rjk, z)/norm_rjk**2
        dt[3*self.j:3*(self.j+1)] = (- dt[3*self.i:3*(self.i+1)]
                                     - dt[3*self.k:3*(self.k+1)])
        return dt


class Dihedral(Primitive):

    def __init__(self, i, j, k, l):  # noqa
        """Implementation follows:
        Blondel, A. and Karplus, M., J. Comput. Chem., 17: 1132-1141. (1996)
        """
        self.i, self.j, self.k, self.l = i, j, k, l  # noqa

    def value(self, xyzs):
        f = xyzs[self.i, :] - xyzs[self.j, :]
        g = xyzs[self.j, :] - xyzs[self.k, :]
        h = xyzs[self.l, :] - xyzs[self.k, :]
        a = np.cross(f, g)
        b = np.cross(h, g)
        norm_g = np.linalg.norm(g)
        w = np.arctan2(
                np.dot(np.cross(b, a), g/norm_g),
                np.dot(a, b))
        return w

    def derivative(self, xyzs):
        """Formula 27 (i,j,k,l)"""
        f = xyzs[self.i, :] - xyzs[self.j, :]
        g = xyzs[self.j, :] - xyzs[self.k, :]
        h = xyzs[self.l, :] - xyzs[self.k, :]
        a = np.cross(f, g)
        b = np.cross(h, g)
        norm_g = np.linalg.norm(g)
        a_sq = np.dot(a, a)
        b_sq = np.dot(b, b)
        dw = np.zeros(xyzs.size)
        dw[3*self.i:3*(self.i+1)] = - norm_g/a_sq*a
        dw[3*self.j:3*(self.j+1)] = (norm_g/a_sq*a
                                     + np.dot(f, g)/(a_sq*norm_g)*a
                                     - np.dot(h, g)/(b_sq*norm_g)*b)
        dw[3*self.k:3*(self.k+1)] = (np.dot(h, g)/(b_sq*norm_g)*b
                                     - np.dot(f, g)/(a_sq*norm_g)*a
                                     - norm_g/b_sq*b)
        dw[3*self.l:3*(self.l+1)] = norm_g/b_sq*b
        return dw


def find_connectivity(atoms, threshold=1.25, connect_fragments=True):
    """Follows the algorithm outlined in Section II A of
    Billeter et al. Phys. Chem. Chem. Phys., 2000, 2, 2177-2186
    """
    bonds = []
    N = len(atoms)
    xyzs = atoms.get_positions()
    types = atoms.get_atomic_numbers()

    r2 = np.zeros((N, N))
    np.fill_diagonal(r2, np.inf)
    neighbors_per_atom = [[] for _ in range(N)]

    for i in range(N):
        for j in range(i+1, N):
            r2[i, j] = r2[j, i] = np.sum((xyzs[i, :] - xyzs[j, :])**2)
            # The paper clearly states that squared distances are compared
            # without squaring the threshold:
            if (r2[i, j] < threshold * (
                    ase.data.covalent_radii[types[i]]
                    + ase.data.covalent_radii[types[j]])**2):
                bonds.append((i, j))
                neighbors_per_atom[i].append(j)
                neighbors_per_atom[j].append(i)

    if connect_fragments:
        # Check for disconnected fragments
        def depth_first_search(fragment, i, visited):
            visited[i] = True
            fragment.append(i)
            for j in neighbors_per_atom[i]:
                if not visited[j]:
                    fragment = depth_first_search(fragment, j, visited)
            return fragment

        visited = [False] * N
        fragments = []
        for i in range(N):
            if not visited[i]:
                fragment = []
                fragments.append(depth_first_search(fragment, i, visited))

        # If there are more than 1 fragment connect shortest distance atoms
        while len(fragments) > 1:
            # Merge first fragment with nearest fragment by adding a bond
            r_min = np.inf
            bond_to_add = ()
            for frag_ind, other_frag in enumerate(fragments[1:]):
                for i in fragments[0]:
                    for j in other_frag:
                        if r2[i, j] < r_min:
                            r_min = r2[i, j]
                            # Plus one because of iterating fragments[1:]
                            fragment_to_merge = frag_ind + 1
                            bond_to_add = tuple(sorted([i, j]))
            bonds.append(bond_to_add)
            fragments[0].extend(fragments.pop(fragment_to_merge))
    return bonds
