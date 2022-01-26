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


def test_xtb_high_spin(tmpdir):
    testName = "xtb_acac_spin5"
    threshMLBL = 0.01
    threshLG = 0.01
    threshOG = 0.01
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_gfnff(tmpdir):
    testName = "gfnff_NH3_BA"
    threshMLBL = 0.01
    threshLG = 0.01
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
