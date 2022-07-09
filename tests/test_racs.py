import pytest
import pickle
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D
from molSimplify.Informatics.RACassemble import create_OHE
from molSimplify.Informatics.lacRACAssemble import get_descriptor_vector
from pkg_resources import resource_filename, Requirement


def test_Fe_CO_6():
    mol = mol3D()
    mol.addAtom(atom3D(Sym='Fe', xyz=[0., 0., 0.]))
    r_FeC = 1.942
    r_CO = 1.125
    for axis in range(3):
        for direction in [1., -1.]:
            xyz_C = [0., 0., 0.]
            xyz_C[axis] = direction * r_FeC
            mol.addAtom(atom3D(Sym='C', xyz=xyz_C))
            xyz_O = [0., 0., 0.]
            xyz_O[axis] = direction * (r_FeC + r_CO)
            mol.addAtom(atom3D(Sym='O', xyz=xyz_O))

    depth = 3
    properties = ['chi', 'Z', 'T', 'S', 'I']
    start_scopes_product = [('f', 'all'), ('mc', 'all'), ('lc', 'ax'),
                            ('lc', 'eq'), ('f', 'ax'), ('f', 'eq')]
    start_scopes_difference = [('mc', 'all'), ('lc', 'ax'), ('lc', 'eq')]
    misc_properties = ['dent', 'charge']
    misc_scopes = ['ax', 'eq']

    names, values = get_descriptor_vector(mol, depth=depth)
    features = {name: val for name, val in zip(names, values)}

    # Product: # start/scopes x # properties x (depth + 1)
    # Difference: # start/scopes X # properties x (depth + 1)
    # Misc: # misc_properties x # misc_scopes
    assert len(features) == (len(start_scopes_product) * len(properties) * (depth + 1)
                             + len(start_scopes_difference) * len(properties) * (depth + 1)
                             + len(misc_properties) * len(misc_scopes))

    # Product RACs
    for start, scope in start_scopes_product:
        for prop in properties:
            for d in range(depth+1):
                assert f'{start}-{prop}-{d}-{scope}' in features

    # Difference RACs. Note the property 'I' and depth 0 are typically not used.
    # See Section 2.2 in J. Phys. Chem. A ,2017, 121, 46, 8939-8954 for details.
    for start, scope in start_scopes_difference:
        for prop in properties:
            for d in range(depth+1):
                assert f'D_{start}-{prop}-{d}-{scope}' in features

    # Misc:
    for prop in misc_properties:
        for scope in misc_scopes:
            assert f'misc-{prop}-{scope}' in features

    assert features == mol.get_features(depth=depth)

    ref_path = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/refs/racs/racs_Fe_(CO)_6.pickle')
    with open(ref_path, 'rb') as fin:
        ref_features = pickle.load(fin)

    assert features.keys() == ref_features.keys()
    for key, val in features.items():
        assert abs(val - ref_features[key]) < 1e-4


def test_Mn_water2_ammonia_furan2_ammonia():
    xyz_path = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/refs/racs/'
        'mn_furan_water_ammonia_furan_water_ammonia.xyz')
    mol = mol3D()
    mol.readfromxyz(xyz_path)
    features = mol.get_features()

    ref_path = resource_filename(
        Requirement.parse('molSimplify'),
        'tests/refs/racs/'
        'racs_Mn_furan_water_ammonia_furan_water_ammonia.pickle')
    with open(ref_path, 'rb') as fin:
        ref_features = pickle.load(fin)

    assert features.keys() == ref_features.keys()
    for key, val in features.items():
        assert abs(val - ref_features[key]) < 1e-4


def test_create_OHE():
    ohe_names, ohe_values = create_OHE('Cr', '2')
    assert ohe_names == ['ox2', 'ox3', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8']
    assert ohe_values == [1, 0, 0, 1, 0, 0, 0, 0]

    _, ohe_values = create_OHE('Cr', '3')
    assert ohe_values == [0, 1, 1, 0, 0, 0, 0, 0]

    _, ohe_values = create_OHE('Mn', '2')
    assert ohe_values == [1, 0, 0, 0, 1, 0, 0, 0]

    _, ohe_values = create_OHE('Mn', '3')
    assert ohe_values == [0, 1, 0, 1, 0, 0, 0, 0]

    _, ohe_values = create_OHE('Fe', '2')
    assert ohe_values == [1, 0, 0, 0, 0, 1, 0, 0]

    _, ohe_values = create_OHE('Fe', '3')
    assert ohe_values == [0, 1, 0, 0, 1, 0, 0, 0]

    _, ohe_values = create_OHE('Co', '2')
    assert ohe_values == [1, 0, 0, 0, 0, 0, 1, 0]

    _, ohe_values = create_OHE('Co', '3')
    assert ohe_values == [0, 1, 0, 0, 0, 1, 0, 0]

    _, ohe_values = create_OHE('Ni', '2')
    assert ohe_values == [1, 0, 0, 0, 0, 0, 0, 1]
