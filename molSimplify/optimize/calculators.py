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
