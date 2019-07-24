import helperFuncs as hp


def test_example_1(tmpdir):
    testName = "one_empty_good"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh, geo_type="one_empty")
    assert passGeo


def test_example_2(tmpdir):
    testName = "one_empty_bad"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh, geo_type="one_empty")
    assert passGeo
