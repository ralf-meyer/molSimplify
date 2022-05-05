from molSimplify.optimize.io import read_molecule, read_terachem_input
from molSimplify.optimize.calculators import TeraChem
from pkg_resources import resource_filename, Requirement


def test_read_molecule():
    in_file = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/optimize/inputs/acac/acac.inp')
    atoms, xyz_file = read_molecule(in_file)
    assert xyz_file.name == 'fe_oct_2_acac_3_s_5_conf_1.xyz'
    assert atoms.get_initial_magnetic_moments().sum() == 4.
    assert atoms.get_initial_charges().sum() == -1.


def test_read_terachem_input():
    in_file = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/optimize/inputs/acac/acac.inp')
    atoms = read_terachem_input(in_file)
    assert atoms.get_initial_magnetic_moments().sum() == 4.
    assert atoms.get_initial_charges().sum() == -1.
    calc = atoms.calc
    assert type(calc) is TeraChem

    params = {'timings': 'yes',
              'maxit': '500',
              'scrdir': './scr',
              'method': 'ub3lyp',
              'basis': 'lacvps_ecp',
              'spinmult': '5',
              'charge': '-1',
              'gpus': '1',
              'scf': 'diis+a',
              'levelshift': 'yes',
              'levelshiftvala': '0.25',
              'levelshiftvalb': '0.25'}
    for key, val in params.items():
        assert calc.parameters[key] == val
