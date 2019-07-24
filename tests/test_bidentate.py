import helperFuncs as hp


def test_bidentate(tmpdir):
    testName = "bidentate"
    threshMLBL = 0.1
    threshLG = 1.0
    threshOG = 2.0
    [passMultiFileCheck, pass_structures] = hp.runtestMulti(
        tmpdir, testName, threshMLBL, threshLG, threshOG)
    assert passMultiFileCheck
    for passStruct in pass_structures:
        [f, passNumAtoms, passMLBL, passLG, passOG, pass_report] = passStruct
        print(f)
        assert passNumAtoms
        assert passMLBL
        assert passLG
        assert passOG
        assert pass_report
