import helperFuncs as hp


def test_xtb_before(tmpdir):
    testName = "xtb_H2O_before"
    threshMLBL = 0.01
    threshLG = 0.01
    threshOG = 2.0
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_xtb_before_after(tmpdir):
    testName = "xtb_imidazole_BA"
    threshMLBL = 0.01
    threshLG = 0.01
    threshOG = 0.05
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_xtb_ANC_fail(tmpdir):
    testName = "xtb_ANC_fail"
    threshMLBL = 0.01
    threshLG = 0.01
    threshOG = 0.01
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_xtb_high_spin(tmpdir):
    testName = "xtb_bipy_spin5"
    threshMLBL = 0.01
    threshLG = 0.05
    threshOG = 0.5
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin


def test_xtb_final_opt(tmpdir):
    testName = "xtb_final_opt"
    threshMLBL = 0.01
    threshLG = 0.01
    threshOG = 0.05
    (passNumAtoms, passMLBL, passLG,
     passOG, pass_report, pass_qcin) = hp.runtest(
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
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
        tmpdir, testName, threshMLBL, threshLG, threshOG, seed=31415)
    assert passNumAtoms
    assert passMLBL
    assert passLG
    assert passOG
    assert pass_report
    assert pass_qcin
