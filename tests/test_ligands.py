import pytest
from molSimplify.Scripts.io import readdict, lig_load
from pkg_resources import resource_filename, Requirement

path = resource_filename(Requirement.parse(
    'molSimplify'), 'molSimplify/Ligands/ligands.dict')
lig_dict = readdict(path)


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
