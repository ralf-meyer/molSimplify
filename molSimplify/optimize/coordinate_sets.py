import numpy as np
from molSimplify.utils.exceptions import ConvergenceError
from molSimplify.optimize.hessian_guess import LindhHessian
from warnings import warn


class CoordinateSet():

    def size(self):
        """Return the number of internal coordinates."""

    def B(self, xyzs):
        """Transformation matrix for the linearized transformation from
        Cartesian to internal coordinates."""

    def to_internals(self, xyzs):
        """Transform a given Cartesian geometry to the internal representation.
        """

    def to_cartesians(self, dq, xyzs_ref):
        """For a given step in internal coordinates find the new Cartesian
        geometry."""

    def diff_internals(self, xyzs1, xyzs2):
        """Calculate the distance (1-2) between the internal representations of
        two Cartesian geometries. This is a separate method as some internal
        representation might need cleaning up, e.g. resticting angle to a
        specific range."""

    def force_to_internals(self, xyzs, force_cart):
        """Transfrom a Cartesian force vector to the internal coordinate
        system."""

    def hessian_to_internals(self, xyzs, hess_cart, grad_cart=None):
        """Transform a Cartesian Hessian matrix to the internal coordinate
        system."""


class CartesianCoordinates(CoordinateSet):

    def __init__(self, atoms):
        self.n_atoms = len(atoms)

    def size(self):
        return 3*self.n_atoms

    def B(self, xyzs):
        return np.eye(xyzs.size)

    def to_internals(self, xyzs):
        return xyzs.flatten()

    def to_cartesians(self, dq, xyzs_ref):
        return xyzs_ref + dq.reshape(xyzs_ref.shape)

    def diff_internals(self, xyzs1, xyzs2):
        return xyzs1.flatten() - xyzs2.flatten()

    def force_to_internals(self, xyzs, force_cart):
        return force_cart.flatten()

    def hessian_to_internals(self, xyzs, hess_cart, grad_cart=None):
        return hess_cart


class InternalCoordinates(CoordinateSet):

    def __init__(self, primitives):
        self.primitives = primitives

    def size(self):
        return len(self.primitives)

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

    def to_cartesians(self, dq, xyzs_ref, tol_q=1e-10, tol_x=1e-10,
                      maxstep=0.05, maxiter=50, recursion_depth=0):
        xyzs = xyzs_ref.copy()
        dq_start = dq.copy()
        step = np.infty * np.ones_like(xyzs)
        for _ in range(maxiter):
            if np.linalg.norm(dq) < tol_q or np.linalg.norm(step) < tol_x:
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
        if recursion_depth >= 3:
            raise ConvergenceError('Transformation to Cartesians not converged'
                                   f' within {maxiter} iterations')
        # Else warn, reduce dq_step by half, and try again
        warn('Reducing step in transformation to Cartesians')
        return self.to_cartesians(dq_start/2, xyzs_ref, tol_q=tol_q,
                                  tol_x=tol_x, maxstep=maxstep,
                                  maxiter=maxiter,
                                  recursion_depth=recursion_depth + 1)

    def diff_internals(self, xyzs1, xyzs2):
        dq = np.zeros(len(self.primitives))
        for i, prim in enumerate(self.primitives):
            dq[i] = prim.diff(xyzs1, xyzs2)
        return dq

    def force_to_internals(self, xyzs, force_cart):
        return self.Binv(xyzs) @ force_cart.flatten()

    def hessian_to_internals(self, xyzs, hess_cart, grad_cart=None):
        Binv = self.Binv(xyzs)

        if grad_cart is not None:
            raise NotImplementedError('Transformation including gradient term '
                                      'is not implemented yet')
            # hess_cart -= Binv @ grad_cart @ self.second_derivatives(xyzs)

        hess_int = Binv @ hess_cart @ Binv.T
        return hess_int


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

    def size(self):
        return self.U.shape[1]

    def B(self, xyzs):
        return self.U.T @ InternalCoordinates.B(self, xyzs)

    def to_internals(self, xyzs):
        return self.U.T @ InternalCoordinates.to_internals(self, xyzs)

    def diff_internals(self, xyzs1, xyzs2):
        return self.U.T @ InternalCoordinates.diff_internals(self, xyzs1,
                                                             xyzs2)


class ApproximateNormalCoordinates(CoordinateSet):

    def __init__(self, atoms, threshold=0.):
        self.threshold = threshold
        self.build(atoms)

    def build(self, atoms):
        H = LindhHessian(h_trans=0., h_rot=0.).build(atoms)
        vals, V = np.linalg.eigh(H)
        self.V = V[:, np.abs(vals) >= self.threshold]
        self.x0 = atoms.get_positions()

    def size(self):
        return self.V.shape[1]

    def B(self, xyzs):
        return self.V.T

    def to_internals(self, xyzs):
        return (xyzs - self.x0).flatten() @ self.V

    def to_cartesians(self, dq, xyzs_ref):
        return xyzs_ref + (self.V @ dq).reshape(xyzs_ref.shape)

    def diff_internals(self, xyzs1, xyzs2):
        return (xyzs1 - xyzs2).flatten() @ self.V

    def force_to_internals(self, xyzs, force_cart):
        return force_cart @ self.V

    def hessian_to_internals(self, xyzs, hess_cart, grad_cart=None):
        return self.V.T @ hess_cart @ self.V
