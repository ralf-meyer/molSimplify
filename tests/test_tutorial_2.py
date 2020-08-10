import helperFuncs as hp


def test_tutorial_2(tmpdir):
    testName = "tutorial_2"
    threshOG = 2.0
    [passNumAtoms, passOG] = hp.runtest_slab(tmpdir, testName, threshOG)
    assert passNumAtoms
    assert passOG