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
        if self.eref is None:
            self._calc_reference(xyzs)
        # Define two directions perpendicular to rik using the reference vector
        u = np.cross(eik, self.eref)
        u /= np.linalg.norm(u)
        # Since eik and u are perpendicular and normalized w is normalized by
        # construction.
        w = np.cross(eik, u)

        rji = xyzs[self.i, :] - xyzs[self.j, :]
        eji = rji / np.linalg.norm(rji)
        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        ejk = rjk / np.linalg.norm(rjk)
        if self.axis == 0:
            return np.dot(eji, u) + np.dot(ejk, u)
        # Else use w as projection axis
        return np.dot(eji, w) + np.dot(ejk, w)

    def derivative(self, xyzs):
        def d_unit_vector(a):
            term1 = np.eye(3)/np.linalg.norm(a)
            term2 = np.outer(a, a)/(np.linalg.norm(a)**3)
            answer = term1-term2
            return answer

        def d_cross_ab(a, b, da, db):
            answer = np.zeros((da.shape[0], 3), dtype=float)
            for i in range(da.shape[0]):
                answer[i] = np.cross(a, db[i]) + np.cross(da[i], b)
            return answer

        dt = np.zeros(xyzs.size)
        rik = xyzs[self.k, :] - xyzs[self.i, :]
        eik = rik / np.linalg.norm(rik)
        if self.eref is None:
            self._calc_reference(xyzs)
        c1 = np.cross(eik, self.eref)
        u = c1 / np.linalg.norm(c1)
        w = np.cross(eik, u)
        rji = xyzs[self.i, :] - xyzs[self.j, :]
        eji = rji / np.linalg.norm(rji)
        rjk = xyzs[self.k, :] - xyzs[self.j, :]
        ejk = rjk / np.linalg.norm(rjk)
        # Derivative terms
        # Derivative of reference vector is zero:
        deref = np.zeros((3, 3))
        deik = d_unit_vector(rik)
        dc1 = d_cross_ab(eik, self.eref, deik, deref)
        du = np.dot(dc1, d_unit_vector(c1))
        dw = d_cross_ab(eik, u, deik, du)
        deji = d_unit_vector(rji)
        dejk = d_unit_vector(rjk)
        if self.axis == 0:
            dt[3*self.i:3*(self.i+1)] = (np.dot(deji, u) + np.dot(-du, eji)
                                         + np.dot(-du, ejk))
            dt[3*self.j:3*(self.j+1)] = np.dot(-deji, u) + np.dot(-dejk, u)
            dt[3*self.k:3*(self.k+1)] = (np.dot(du, eji) + np.dot(du, ejk)
                                         + np.dot(dejk, u))
        else:
            dt[3*self.i:3*(self.i+1)] = (np.dot(deji, w) + np.dot(-dw, eji)
                                         + np.dot(-dw, ejk))
            dt[3*self.j:3*(self.j+1)] = np.dot(-deji, w) + np.dot(-dejk, w)
            dt[3*self.k:3*(self.k+1)] = (np.dot(dw, eji) + np.dot(dw, ejk)
                                         + np.dot(dejk, w))
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
