import numpy as np
from molSimplify.utils.exceptions import ConvergenceError


class CartesianCoordinates():

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


class InternalCoordinates():

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
                      maxstep=0.05, maxiter=50):
        xyzs = xyzs_ref.copy()
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
        raise ConvergenceError('Transformation to Cartesians not converged '
                               f'within {maxiter} iterations')

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