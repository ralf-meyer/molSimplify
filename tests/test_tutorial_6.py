import helperFuncs as hp


def test_tutorial_6(tmpdir):
    testName = "tutorial_6"
    threshOG = 2.0
    [passNumAtoms, passOG] = hp.runtest_molecule_on_slab(tmpdir, testName, threshOG)
    assert passNumAtoms
    assert passOG
