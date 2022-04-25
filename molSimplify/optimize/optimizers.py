import numpy as np
import ase.optimize


class BFGS(ase.optimize.BFGS):
    """Adaptation of ASEs implementation of the BFGS optimizer to allow for
    arbitrary (internal) coordinate systems.
    """

    def __init__(self, atoms, coordinate_set, maxstep_internal=1.0, **kwargs):
        if kwargs.get('use_line_search', False):
            raise NotImplementedError('Line search is not implemented yet.')
        self.coord_set = coordinate_set
        ase.optimize.BFGS.__init__(self, atoms, **kwargs)
        self.maxstep_internal = maxstep_internal

    def initialize(self):
        # initial hessian
        # MOD: initialize to correct size
        self.H0 = np.eye(self.coord_set.size()) * self.alpha

        self.H = None
        self.r0 = None
        self.f0 = None

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        r = atoms.get_positions()
        # MOD: Transform forces to internal coordinates
        f = self.coord_set.force_to_internals(r, f.flatten())
        # MOD: remove the flattening here
        self.update(r, f, self.r0, self.f0)
        omega, V = np.linalg.eigh(self.H)

        # FUTURE: Log this properly
        # # check for negative eigenvalues of the hessian
        # if any(omega < 0):
        #     n_negative = len(omega[omega < 0])
        #     msg = '\n** BFGS Hessian has {} negative eigenvalues.'.format(
        #         n_negative
        #     )
        #     print(msg, flush=True)
        #     if self.logfile is not None:
        #         self.logfile.write(msg)
        #         self.logfile.flush()

        # MOD: step in internals
        dq = np.dot(V, np.dot(f, V) / np.fabs(omega))
        if np.max(np.abs(dq)) > self.maxstep_internal:
            dq *= self.maxstep_internal / np.max(np.abs(dq))
        # Transform to Cartesians
        dr = self.coord_set.to_cartesians(dq, r) - r
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)
        self.r0 = r.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))

    def update(self, r, f, r0, f0):
        if self.H is None:
            self.H = self.H0
            return
        dr = self.coord_set.diff_internals(r, r0).flatten()

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        df = f - f0
        a = np.dot(dr, df)
        if a > 0:
            print('Skipping BFGS update to conserve positive '
                  'definiteness.')
            return
        dg = np.dot(self.H, dr)
        b = np.dot(dr, dg)
        self.H -= np.outer(df, df) / a + np.outer(dg, dg) / b


class RFO(BFGS):

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        r = atoms.get_positions()
        f = self.coord_set.force_to_internals(r, f.flatten())
        self.update(r, f, self.r0, self.f0)

        # extended Hessian matrix
        H_ext = np.block([[self.H, -f[:, np.newaxis]], [-f, 0.]])

        _, V = np.linalg.eigh(H_ext)

        # Step is calculated by proper rescaling of the eigenvector
        # corresponding to the lowest (first) eigenvalue.
        dq = V[:-1, 0] / V[-1, 0]
        if np.max(np.abs(dq)) > self.maxstep_internal:
            dq *= self.maxstep_internal / np.max(np.abs(dq))
        # Transform to Cartesians
        dr = self.coord_set.to_cartesians(dq, r) - r
        steplengths = (dr**2).sum(1)**0.5
        dr = self.determine_step(dr, steplengths)
        atoms.set_positions(r + dr)
        self.r0 = r.copy()
        self.f0 = f.copy()
        self.dump((self.H, self.r0, self.f0, self.maxstep))


class LBFGS(ase.optimize.LBFGS):
    """Adaptation of ASEs implementation of the LBFGS optimizer to allow for
    arbitrary (internal) coordinate systems.
    """

    def __init__(self, atoms, coordinate_set, maxstep_internal=1.0, **kwargs):
        if kwargs.get('use_line_search', False):
            raise NotImplementedError('Line search is not implemented yet.')
        ase.optimize.LBFGS.__init__(self, atoms, **kwargs)
        self.coord_set = coordinate_set
        self.maxstep_internal = maxstep_internal

    def step(self, f=None):
        """Take a single step

        This is where actual modifications to account for internal coordinates
        have to be made. Modifications are highlighted by a # MOD comment.

        Use the given forces, update the history and calculate the next step --
        then take it"""

        if f is None:
            f = self.atoms.get_forces()

        r = self.atoms.get_positions()

        # MOD: Transform forces to internal coordinates
        f = self.coord_set.force_to_internals(r, f.flatten())
        # print('q internal', self.coord_set.to_internals(r))
        # print('f internal', f)

        self.update(r, f, self.r0, self.f0)

        s = self.s
        y = self.y
        rho = self.rho
        H0 = self.H0

        loopmax = np.min([self.memory, self.iteration])
        a = np.empty((loopmax,), dtype=np.float64)

        # ## The algorithm itself:
        q = -f.reshape(-1)
        for i in range(loopmax - 1, -1, -1):
            a[i] = rho[i] * np.dot(s[i], q)
            q -= a[i] * y[i]
        z = H0 * q

        for i in range(loopmax):
            b = rho[i] * np.dot(y[i], z)
            z += s[i] * (a[i] - b)

        # MOD: Calculate internal step
        dq = -z
        if np.max(np.abs(dq)) > self.maxstep_internal:
            dq *= self.maxstep_internal / np.max(np.abs(dq))
        # MOD: Transform internal step to Cartesian step
        self.p = self.coord_set.to_cartesians(dq, r) - r
        # ##

        g = -f
        if self.use_line_search is True:
            e = self.func(r)
            self.line_search(r, g, e)
            dr = (self.alpha_k * self.p).reshape(len(self.atoms), -1)
        else:
            self.force_calls += 1
            self.function_calls += 1
            dr = self.determine_step(self.p) * self.damping
        # MOD: xyzs instead of r
        self.atoms.set_positions(r + dr)

        self.iteration += 1
        self.r0 = r
        self.f0 = -g
        self.dump((self.iteration, self.s, self.y,
                   self.rho, self.r0, self.f0, self.e0, self.task))

    def update(self, r, f, r0, f0):
        """Update everything that is kept in memory

        This function is mostly here to allow for replay_trajectory.
        """
        if self.iteration > 0:
            # MOD: step s should be calculated in internals
            s0 = self.coord_set.diff_internals(r, r0)
            self.s.append(s0)

            # We use the gradient which is minus the force!
            y0 = f0.reshape(-1) - f.reshape(-1)
            self.y.append(y0)

            rho0 = 1.0 / np.dot(y0, s0)
            self.rho.append(rho0)

        if self.iteration > self.memory:
            self.s.pop(0)
            self.y.pop(0)
            self.rho.pop(0)

    def replay_trajectory(self, traj):
        raise NotImplementedError('Trajectory replay is not yet supported!')
