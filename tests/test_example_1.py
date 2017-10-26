import helperFuncs as hp

def test_example_1(tmpdir):
    [pass_xyz,pass_report] = hp.runtest(tmpdir,"example_1",0.1)
    assert pass_xyz and pass_report
