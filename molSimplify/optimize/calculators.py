import io
import numpy as np
import ase.io
import ase.units
import ase.calculators.calculator as ase_calculator
try:  # For compatibility with openbabel < 3.0
    from openbabel import openbabel
except ImportError:
    import openbabel


_xtb_methods = ['xtb']
_openbabel_methods = ['uff', 'mmff94', 'gaff']
_available_methods = _xtb_methods + _openbabel_methods


def get_calculator(method):
    if method.lower() in _xtb_methods:
        from xtb.ase.calculator import XTB
        methods = dict(xtb=('GFN2-xTB', 300.0),
                       gfnff=('GFNFF', 0.0))
        # TODO: can't get GFNFF to work!
        return XTB(method=methods[method][0],
                   electronic_temperature=methods[method][1])
    elif method.lower() in _openbabel_methods:
        return OpenbabelFF(ff=method.upper())


class OpenbabelFF(ase_calculator.Calculator):

    implemented_properties = ['energy', 'forces']

    nolabel = True

    default_parameters = {'ff': 'UFF'}

    ob_units = {'kcal/mol': ase.units.kcal/ase.units.mol,
                'kJ/mol': ase.units.kJ/ase.units.mol}

    def __init__(self, **kwargs):
        ase_calculator.Calculator.__init__(self, **kwargs)
        self.outputname = 'openbabelff'

    def initialize(self, atoms):
        obMol = OpenbabelFF.atoms2OBmol(atoms)
        self.ff = openbabel.OBForceField.FindForceField(self.parameters.ff)
        if not self.ff.Setup(obMol):
            for atom in openbabel.OBMolAtomIter(obMol):
                print(atom.GetAtomicNum(), atom.GetFormalCharge())
            raise RuntimeError(
                f'Could not setup force field {self.parameters.ff}')
        self.energy_unit = self.ob_units[self.ff.GetUnit()]

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=ase_calculator.all_changes):
        ase_calculator.Calculator.calculate(self, atoms, properties,
                                            system_changes)

        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        obMol = OpenbabelFF.atoms2OBmol(self.atoms)
        self.ff.SetCoordinates(obMol)
        self.results['energy'] = self.ff.Energy()*self.energy_unit
        grad = np.zeros((len(atoms), 3))
        for i, atom_i in enumerate(openbabel.OBMolAtomIter(obMol)):
            grad_i = self.ff.GetGradient(atom_i)
            grad[i, :] = (grad_i.GetX(), grad_i.GetY(), grad_i.GetZ())
        # Note, GetGradient() returns the negative gradient, so no sign
        # inversion needed here.
        self.results['forces'] = grad * self.energy_unit/ase.units.Ang

    @staticmethod
    def atoms2OBmol(atoms):
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat('xyz')
        obMol = openbabel.OBMol()

        with io.StringIO() as stream:
            ase.io.write(stream, atoms, format='xyz')
            obConversion.ReadString(obMol, stream.getvalue())
        return obMol


class TwoDCalculator(ase_calculator.Calculator):
    """Base class for two dimensional benchmark systems."""

    implemented_properties = ['energy', 'forces']

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=ase_calculator.all_changes):
        ase_calculator.Calculator.calculate(self, atoms, properties,
                                            system_changes)

        xyzs = self.atoms.get_positions()
        x = xyzs[0, 0]
        dx = np.zeros_like(xyzs)
        dx[0, 0] = 1.0
        y = xyzs[0, 1]
        dy = np.zeros_like(xyzs)
        dy[0, 1] = 1.0

        self.results['energy'] = self.energy(x, y)
        gx, gy = self.gradient(x, y)
        self.results['forces'] = -gx * dx - gy * dy


class CerjanMillerSurface(TwoDCalculator):
    a = 1
    b = 1
    c = 1

    def energy(self, x, y):
        return ((self.a - self.b * y**2) * x**2 * np.exp(-x**2)
                + (self.c/2) * y**2)

    def gradient(self, x, y):
        exp = np.exp(-x**2)
        return (2 * (self.a - self.b * y**2) * exp * x * (1 - x**2),
                - 2 * self.b * y * exp * x**2 + self.c * y)


class AdamsSurface(TwoDCalculator):

    def energy(self, x, y):
        return (2 * x**2 * (4 - x) + y**2 * (4 + y)
                - x*y * (6 - 17 * np.exp(-(x**2 + y**2)/4)))

    def gradient(self, x, y):
        dx1 = (-17/2 * (x**2 - 2) * y * np.exp(-(x**2 + y**2)/4)
               - 6 * x**2 + 16 * x - 6 * y)
        dx2 = (-17/2 * (y**2 - 2) * x * np.exp(-(x**2 + y**2)/4)
               + 3 * y**2 + 8 * y - 6 * x)
        return dx1, dx2


class MuellerBrownSurface(TwoDCalculator):
    A = np.array([-200, -100, -170, 15])
    a = np.array([-1, -1, -6.5, 0.7])
    b = np.array([0, 0, 11, 0.6])
    c = np.array([-10, -10, -6.5, 0.7])
    x0 = np.array([1, 0, -0.5, -1])
    y0 = np.array([0, 0.5, 1.5, 1])

    def _v(self, x, y):
        return self.A * np.exp(self.a*(x - self.x0)**2
                               + self.b*(x - self.x0)*(y - self.y0)
                               + self.c*(y - self.y0)**2)

    def energy(self, x, y):
        return np.sum(self._v(x, y))

    def gradient(self, x, y):
        v = self._v(x, y)
        return ((2*self.a*(x - self.x0) + self.b*(y - self.y0)) @ v,
                (self.b*(x - self.x0) + 2*self.c*(y - self.y0)) @ v)
