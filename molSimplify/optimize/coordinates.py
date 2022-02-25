import numpy as np
from abc import abstractmethod


class Primitive():

    @abstractmethod
    def value(self, xyzs):
        ...

    @abstractmethod
    def derivative(self, xyzs):
        ...

    def diff(self, xyzs1, xyzs2):
        return self.value(xyzs1) - self.value(xyzs2)


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
        return 1./Distance.value(self, xyzs)

    def derivative(self, xyzs):
        rij = xyzs[self.i, :] - xyzs[self.j, :]
        q = 1.0/np.linalg.norm(rij)
        dq = np.zeros(xyzs.size)
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


class LinearAngle(Primitive):
    def __init__(self, i, j, k, axis):
        """Closely follows the implementation in geomeTRIC

        Parameters
        ----------
        i : int
            Index of the first atom
        j : int
            Index of the center atom
        k : int
            Index of the last atom
        axis : int
            Projection axis. Can take values 0 or 1.
        """
        self.i, self.j, self.k = i, j, k
        self.axis = axis
        self.eref = None

    def _calc_reference(self, xyzs):
        rik = xyzs[self.k, :] - xyzs[self.i, :]
        # Cartesian axes.
        cart_vecs = np.eye(3)
        # Select Cartesian axis with the least overlap with rik as
        # reference direction.
        ind = np.argmin([np.dot(ei, rik)**2 for ei in cart_vecs])
        self.eref = cart_vecs[ind]

    def value(self, xyzs):
        # Unit vector pointing from i to k.
        rik = xyzs[self.k, :] - xyzs[self.i, :]
        eik = rik / np.linalg.norm(rik)
        rji = xyzs[self.i, :] - xyzs[self.j, :]
        eji = rji / np.linalg.norm(rji)
        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        ejk = rjk / np.linalg.norm(rjk)

        if self.eref is None:
            self._calc_reference(xyzs)
        # Define the vector u perpendicular to rik using the reference vector
        u = np.cross(eik, self.eref)
        u /= np.linalg.norm(u)

        if self.axis == 0:
            return np.dot(eji, u) + np.dot(ejk, u)
        # Else use a vector w perpendicular to rik and u as projection axis.
        # Since eik and u are perpendicular and normalized w is normalized by
        # construction.
        w = np.cross(eik, u)
        return np.dot(eji, w) + np.dot(ejk, w)

    def derivative(self, xyzs):
        # Initialize return array
        dt = np.zeros(xyzs.size)

        rik = xyzs[self.k, :] - xyzs[self.i, :]
        norm_ik = np.linalg.norm(rik)
        eik = rik / norm_ik
        # deik/drik
        deik = np.eye(3) / norm_ik - np.outer(rik, rik) / norm_ik**3

        rji = xyzs[self.i, :] - xyzs[self.j, :]
        norm_ji = np.linalg.norm(rji)
        eji = rji / norm_ji
        # deji/drji
        deji = np.eye(3) / norm_ji - np.outer(rji, rji) / norm_ji**3

        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        norm_jk = np.linalg.norm(rjk)
        ejk = rjk / norm_jk
        # dejk/drjk
        dejk = np.eye(3) / norm_jk - np.outer(rjk, rjk) / norm_jk**3

        # Setup first projection axis u
        if self.eref is None:
            self._calc_reference(xyzs)
        u_raw = np.cross(eik, self.eref)
        # Since eref is constant: deref/drik = 0. Caution: While other
        # derivative matrices defined here are symmetric du_raw is not!
        du_raw = np.cross(deik, self.eref, axis=0)
        # Normalization
        norm_u = np.linalg.norm(u_raw)
        u = u_raw / norm_u
        # Inner derivative of norm_u in the second term is again du_raw
        du = du_raw / norm_u - np.outer(u_raw, u_raw) @ du_raw / norm_u**3

        if self.axis == 0:
            # derivative w.r.t. atom i: drji/dri = 1, drik/dri = -1
            dt[3*self.i:3*(self.i+1)] = (np.dot(deji, u) + np.dot(eji, -du)
                                         + np.dot(ejk, -du))
            # derivative w.r.t. atom j: drji/drj = -1, drjk/drj = -1
            # u is independent of atom j : du/drj = 0
            dt[3*self.j:3*(self.j+1)] = np.dot(-deji, u) + np.dot(-dejk, u)
            # derivative w.r.t atom k: drik/drk = 1, drjk/drk = 1
            dt[3*self.k:3*(self.k+1)] = (np.dot(dejk, u) + np.dot(eji, du)
                                         + np.dot(ejk, du))
        else:
            # Setup second projection axis
            w = np.cross(eik, u)
            # Derivative w.r.t rik
            dw = np.cross(deik, u, axis=0) + np.cross(eik, du, axis=0)
            # derivative w.r.t. atom i: drji/dri = 1, drik/dri = -1
            dt[3*self.i:3*(self.i+1)] = (np.dot(deji, w) + np.dot(eji, -dw)
                                         + np.dot(ejk, -dw))
            # derivative w.r.t. atom j: drji/drj = -1, drjk/drj = -1
            # w is independent of the atom j: dw/drj = 0
            dt[3*self.j:3*(self.j+1)] = np.dot(-deji, w) + np.dot(-dejk, w)
            # derivative w.r.t atom k: drik/drk = 1, drjk/drk = 1
            dt[3*self.k:3*(self.k+1)] = (np.dot(dejk, w) + np.dot(eji, dw)
                                         + np.dot(ejk, dw))
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

    def diff(self, xyzs1, xyzs2):
        dw = self.value(xyzs1) - self.value(xyzs2)
        if dw > np.pi:
            return dw - 2*np.pi
        if dw < -np.pi:
            return dw + 2*np.pi
        return dw


class Improper(Dihedral):
    """Alias for Dihedral since it is often necessary to distinguish between
    actual dihedrals and improper (out-of-plane) bends."""


class Octahedral(Primitive):

    def __init__(self, a0, a1, a2, a3, a4, a5, a6):
        """ Geometry:
               a5  a4
                | /
                |/
        a1 - - a0 - - a3
               /|
              / |
            a2  a6
        """
        self.a0, self.a1, self.a2, self.a3, self.a4, self.a5, self.a6 = (
            a0, a1, a2, a3, a4, a5, a6)


class InternalCoordinates():

    def __init__(self, primitives):
        self.primitives = primitives

    def B(self, xyzs):
        B = np.zeros((len(self.primitives), xyzs.size))
        for i, prim in enumerate(self.primitives):
            B[i, :] = prim.derivative(xyzs)
        return B

    def Ginv(self, xyzs):
        B = self.B(xyzs)
        return np.linalg.pinv(B @ B.T)

    def Binv(self, xyzs):
        B = self.B(xyzs)
        return np.linalg.pinv(B @ B.T) @ B

    def to_internals(self, xyzs):
        q = np.zeros(len(self.primitives))
        for i, prim in enumerate(self.primitives):
            q[i] = prim.value(xyzs)
        return q

    def to_cartesians(self, dq, xyzs_ref, tol=1e-10, maxstep=0.1, maxiter=50):
        xyzs = xyzs_ref.copy()
        for _ in range(maxiter):
            if np.linalg.norm(dq) < tol:
                return xyzs
            step = (self.Binv(xyzs).T @ dq).reshape(xyzs.shape)
            steplengths = np.sqrt(np.sum(step**2, axis=-1))
            maxlength = np.max(steplengths)
            if maxlength > maxstep:
                step *= maxstep / maxlength
            # Calculate what this step corresponds to in internals
            step_q = self.diff_internals(xyzs + step, xyzs)
            xyzs = xyzs + step
            # Calculate the step for the next iteration
            dq -= step_q
        raise RuntimeError('Transformation to Cartesians not converged within '
                           f'{maxiter} iterations')

    def diff_internals(self, xyzs1, xyzs2):
        dq = np.zeros(len(self.primitives))
        for i, prim in enumerate(self.primitives):
            dq[i] = prim.diff(xyzs1, xyzs2)
        return dq


class DelocalizedCoordinates(InternalCoordinates):

    def __init__(self, primitives, xyzs, threshold=1e-10):
        InternalCoordinates.__init__(self, primitives)
        self.threshold = threshold
        self.delocalize(xyzs)

    def delocalize(self, xyzs):
        B = InternalCoordinates.B(self, xyzs)
        # G matrix without mass weighting
        G = B @ B.T
        w, v = np.linalg.eigh(G)
        # Set of nonredundant eigenvectors (eigenvalue =/= 0)
        self.U = v[:, np.abs(w) > self.threshold]

    def B(self, xyzs):
        return self.U.T @ InternalCoordinates.B(self, xyzs)

    def to_internals(self, xyzs):
        return self.U.T @ InternalCoordinates.to_internals(self, xyzs)

    def diff_internals(self, xyzs1, xyzs2):
        return self.U.T @ InternalCoordinates.diff_internals(self, xyzs1,
                                                             xyzs2)
