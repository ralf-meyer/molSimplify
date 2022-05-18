import pytest
from molSimplify.Scripts.io import readdict, lig_load
from pkg_resources import resource_filename, Requirement
from os import listdir
from os.path import isfile, join

path_dict = resource_filename(Requirement.parse(
    'molSimplify'), 'molSimplify/Ligands/ligands.dict')
lig_dict = readdict(path_dict)


@pytest.mark.parametrize('lig_name', lig_dict.keys())
def test_ligands_dict(lig_name):
    # Assert that the dict entry has all 6 fields
    assert len(lig_dict[lig_name]) == 6
    lig, emsg = lig_load(lig_name)
    # Assert that the ligand could be loaded
    assert emsg is False
    # Assert that the charge of the loaded ligand equals the
    # charge noted in ligands.dict
    charge = int(lig_dict[lig_name][-1][0])
    assert lig.charge == charge
    connecting_atoms = lig_dict[lig_name][-4]
    group = lig_dict[lig_name][-3]
    if 'bidentate' in group:
        assert len(connecting_atoms) == 2
    elif 'tridentate' in group:
        assert len(connecting_atoms) == 3
    elif 'tetradentate' in group:
        assert len(connecting_atoms) == 4


def test_no_repeats():
    # This test ensures no key is used more than once in ligands.dict.
    listed_keys = []
    with open(path_dict) as f:
        array_content = f.readlines()
    array_content = [line for line in array_content if not line.startswith('#')] # Ignoring comments.
    listed_keys = [line.split(':')[0] for line in array_content] # Grabbing the ligand names.
    unique_listed_keys = [item for idx, item in enumerate(listed_keys) if item not in listed_keys[:idx]] # The unique ligand names.
    assert len(listed_keys) == len(set(listed_keys)) # The set will eliminate repeats.


path_folder = resource_filename(Requirement.parse(
    'molSimplify'), 'molSimplify/Ligands/')
# Identify the files in the Ligands folder.
my_files = [f for f in listdir(path_folder) if isfile(join(path_folder, f))]
# Keep only the files that end in .mol or .xyz
my_files = [f for f in my_files if f.endswith('.xyz') or f.endswith('.mol')]
# Getting the geometries listed in ligands.dict
geometries = [entry[1][0] for entry in list(lig_dict.items())]# Requires some ugly indexing to get this done.
@pytest.mark.parametrize('file_name', my_files)
def test_unused_geometries(file_name):
    # This test checks whether any geometry files in the Ligands folder are not used in ligands.dict.
    # This indicates that the geometry file should be removed, as it is clutter.
    assert file_name in geometries


