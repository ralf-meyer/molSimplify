import helperFuncs as hp


def test_tutorial_qm9_part_one(tmpdir):
    testName = "tutorial_qm9_part_one"
    threshOG = 2.0
    [passNumAtoms, passOG] = hp.runtest_num_atoms_in_xyz(xyzfile)
    assert passNumAtoms
    assert passOG
