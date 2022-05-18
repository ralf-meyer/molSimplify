import pytest
import helperFuncs as hp


@pytest.mark.skip
def test_tutorial_qm9_part_one(tmpdir):
    testName = "tutorial_qm9_part_one"
    hp.runtest_num_atoms_in_xyz(tmpdir, testName)

