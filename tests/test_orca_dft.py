import helperFuncs as hp


def test_orca_dft(tmpdir):
    testName = "orca_dft"
    threshMLBL = 0.1
    threshLG = 1.0
    threshOG = 8.0
    [passNumAtoms, passMLBL, passLG, passOG, pass_report, pass_qcin] = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
