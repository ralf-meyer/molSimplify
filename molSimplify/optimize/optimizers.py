from abc import abstractmethod
import time
import numpy as np
import ase.optimize
import ase.units


class ConvergenceMixin():
    """Mixin class used to replace the convergence criteria of an
    molSimplify.optimize.Optimizer class. Usage:
    >>> class NewOptimizer(ConvergenceMixin, OldOptimizer):
            pass
    >>> opt = NewOptimizer(atoms)
    Note: Because of Pythons method resolution order the mixin class
    needs to come first!
    """
    threshold_energy = 1e-4
    threshold_max_step = 1e-3
    threshold_rms_step = 5e-2
    threshold_max_grad = 1e-3
    threshold_rms_grad = 1e-3

    def convergence_condition(self, energy_change, max_step, rms_step,
                              max_grad, rms_grad):
        """This is separated out as some convergence criteria might
        want to implement a more sophisticated logic than just all True"""
        conditions = (energy_change < self.threshold_energy,
                      max_step < self.threshold_max_step,
                      rms_step < self.threshold_rms_step,
                      max_grad < self.threshold_max_grad,
                      rms_grad < self.threshold_rms_grad)
        return all(conditions), conditions

    def irun(self, steps=None):
        """ remove fmax from the argument list"""
        if steps:
            self.max_steps = steps
        return ase.optimize.optimize.Dynamics.irun(self)

    def run(self, steps=None):
        """ remove fmax from the argument list"""
        if steps:
            self.max_steps = steps
        return ase.optimize.optimize.Dynamics.run(self)

    def converged(self, forces=None):
        """Did the optimization converge?"""
        if forces is None:
            forces = self.atoms.get_forces()

        energy = self.atoms.get_potential_energy()
        previous_energy = getattr(self, 'previous_energy', energy)
        energy_change = abs(energy - previous_energy)

        positions = self.atoms.get_positions()
        previous_positions = getattr(self, 'previous_positions',
                                     positions)
        step = positions - previous_positions
        max_step = np.max(np.abs(step))
        rms_step = np.sqrt(np.mean(step**2))

        max_grad = np.max(np.abs(forces))
        rms_grad = np.sqrt(np.mean(forces**2))

        self.previous_energy = energy
        self.previous_positions = positions
        return self.convergence_condition(energy_change, max_step,
                                          rms_step, max_grad, rms_grad)[0]

    def log(self, forces=None):
        if forces is None:
            forces = self.atoms.get_forces()
        max_grad = np.max(np.abs(forces))
        rms_grad = np.sqrt(np.mean(forces**2))
        e = self.atoms.get_potential_energy(
            force_consistent=self.force_consistent)
        e_old = getattr(self, 'previous_energy', e)
        delta_e = e - e_old

        positions = self.atoms.get_positions()
        previous_positions = getattr(self, 'previous_positions', positions)
        step = positions - previous_positions
        max_step = np.max(np.abs(step))
        rms_step = np.sqrt(np.mean(step**2))

        conditions = self.convergence_condition(
            abs(delta_e), max_step, rms_step, max_grad, rms_grad)[1]
        conv = ['*' if c else ' ' for c in conditions]

        T = time.localtime()
        if self.logfile is not None:
            name = self.__class__.__name__
            if self.nsteps == 0:
                args = (" " * len(name), "Step", "Time", "Energy", "delta_e",
                        "grad_max", "grad_rms", "step_max", "step_rms")
                msg = "%s  %4s %8s %15s %15s  %15s  %15s  %15s  %15s\n" % args
                self.logfile.write(msg)

            msg = (f'{name}:  {self.nsteps:3d} {T[3]:02d}:{T[4]:02d}'
                   f':{T[4]:02d} {e:15.6f} {delta_e:15.6f}{conv[0]} '
                   f'{max_grad:15.6f}{conv[1]} {rms_grad:15.6f}{conv[2]} '
                   f'{max_step:15.6f}{conv[3]} {rms_step:15.6f}{conv[4]}\n')
            self.logfile.write(msg)

            self.logfile.flush()


class TerachemConvergence(ConvergenceMixin):
    threshold_energy = 1e-6 * ase.units.Hartree
    threshold_max_step = 1.8e-3 * ase.units.Bohr
    threshold_rms_step = 1.2e-3 * ase.units.Bohr
    threshold_max_grad = 4.5e-4 * ase.units.Hartree / ase.units.Bohr
    threshold_rms_grad = 3.0e-4 * ase.units.Hartree / ase.units.Bohr


class HessianApproximation(np.ndarray):

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj

    @abstractmethod
    def update(self, dr, dg):
        """Update the Hessian using the step dr and change in gradient dg"""


class BFGSHessian(HessianApproximation):

    def update(self, dr, dg):
        a = np.dot(dr, dg)
        if a < 0:
            print('Skipping BFGS update to conserve positive '
                  'definiteness.')
            return
        v = np.dot(self, dr)
        b = np.dot(dr, v)
        self += np.outer(dg, dg) / a - np.outer(v, v) / b


class InternalCoordinatesOptimizer(ase.optimize.optimize.Optimizer):

    defaults = {**ase.optimize.optimize.Optimizer.defaults,
                'maxstep_internal': 1.0, 'H0': 70.0}

    def __init__(self, atoms, coordinate_set, restart=None, logfile='-',
                 trajectory=None, master=None, H0=None,
                 maxstep=None, maxstep_internal=None):

        self.coord_set = coordinate_set
        if H0 is None:
            self.H0 = self.defaults['H0']
        else:
            self.H0 = H0

        if maxstep is None:
            self.maxstep = self.defaults['maxstep']
        else:
            self.maxstep = maxstep
        if maxstep_internal is None:
            self.maxstep_internal = self.defaults['maxstep_internal']
        else:
            self.maxstep_internal = maxstep_internal
        ase.optimize.optimize.Optimizer.__init__(self, atoms, restart, logfile,
                                                 trajectory, master)

    def initialize(self):
        if np.size(self.H0) == 1:
            H0 = np.eye(self.coord_set.size()) * self.H0
        else:
            H0 = self.H0
        self.H = self.hessian_approx(H0)
        self.r0 = None
        self.f0 = None

    def step(self, f=None):
        if f is None:
            f = self.atoms.get_forces()

        r = self.atoms.get_positions()
        # MOD: Transform forces to internal coordinates
        f = self.coord_set.force_to_internals(r, f.flatten())
        # MOD: remove the flattening here
        self.update(r, f, self.r0, self.f0)

        # MOD: step in internals
        dq = self.internal_step(f)
        # Rescale
        maxsteplength_internal = np.max(np.abs(dq))
        if maxsteplength_internal > self.maxstep_internal:
            dq *= self.maxstep_internal / maxsteplength_internal

        # Transform to Cartesians
        dr = self.coord_set.to_cartesians(dq, r) - r
        # Rescale
        maxsteplength = np.max(np.sqrt(np.sum(dr**2, axis=1)))
        if maxsteplength > self.maxstep:
            dr *= self.maxstep / maxsteplength

        self.atoms.set_positions(r + dr)
        self.r0 = r.copy()
        self.f0 = f.copy()
        self.dump((self.coord_set, self.H, self.r0, self.f0, self.maxstep,
                   self.maxstep_internal))

    @abstractmethod
    def internal_step(self, f):
        """this needs to be implemented by subclasses"""

    def read(self):
        (self.coord_set, self.H, self.r0, self.f0, self.maxstep,
         self.maxstep_internal) = self.load()

    def update(self, r, f, r0, f0):
        if r0 is None or f0 is None:  # No update on the first iteration
            return

        dr = self.coord_set.diff_internals(r, r0).flatten()

        if np.abs(dr).max() < 1e-7:
            # Same configuration again (maybe a restart):
            return

        dg = -(f - f0)
        self.H.update(dr, dg)


class BFGS(InternalCoordinatesOptimizer):
    hessian_approx = BFGSHessian

    def internal_step(self, f):
        omega, V = np.linalg.eigh(self.H)
        print(omega)
        return np.dot(V, np.dot(f, V) / np.fabs(omega))


class RFO(InternalCoordinatesOptimizer):

    def __init__(self, *args, mu=0, **kwargs):
        self.mu = mu
        if self.mu == 0:
            self.hessian_approx = BFGSHessian
        else:
            raise NotImplementedError()
            # self.hessian_approx = BofillHessian
        InternalCoordinatesOptimizer.__init__(self, *args, **kwargs)

    def internal_step(self, f):
        # extended Hessian matrix
        H_ext = np.block([[self.H, -f[:, np.newaxis]], [-f, 0.]])

        _, V = np.linalg.eigh(H_ext)

        # Step is calculated by proper rescaling of the eigenvector
        # corresponding to the mu-th eigenvalue. For minimizations
        # the lowest i.e. zeroth eigenvalue is chosen.
        return V[:-1, self.mu] / V[-1, self.mu]


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
        # MOD: Removed option for linesearch
        dr = self.determine_step(self.p) * self.damping
        # MOD: xyzs instead of r
        self.atoms.set_positions(r + dr)

        self.iteration += 1
        self.force_calls += 1
        self.function_calls += 1
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
