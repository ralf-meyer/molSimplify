from molSimplify.optimize.main import read_molecule
from pkg_resources import resource_filename, Requirement


def test_read_molecule():
    in_file = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/optimize/inputs/acac/acac.inp')
    atoms, xyz_file = read_molecule(in_file)
    assert xyz_file.name == 'fe_oct_2_acac_3_s_5_conf_1.xyz'
    assert atoms.get_initial_magnetic_moments().sum() == 4.
    assert atoms.get_initial_charges().sum() == -1.
