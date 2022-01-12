import helperFuncs as hp


def test_xtb_before(tmpdir):
    testName = "xtb_H2O_before"
    threshMLBL = 0.01
    threshLG = 1e-3
    threshOG = 2.0
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_xtb_before_after(tmpdir):
    testName = "xtb_CO_before_after"
    threshMLBL = 0.01
    threshLG = 1e-3
    threshOG = 1e-3
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin
