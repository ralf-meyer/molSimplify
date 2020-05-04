import helperFuncs as hp


def test_example_1(tmpdir):
    testName = "all_flying_away"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_2(tmpdir):
    testName = "broken_ligands"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_3(tmpdir):
    testName = "catom_change"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_4(tmpdir):
    testName = "H_transfer"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_5(tmpdir):
    testName = "ligand_assemble"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_6(tmpdir):
    testName = "ligand_bent"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_7(tmpdir):
    testName = "linear_broken"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_8(tmpdir):
    testName = "methane_trans"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_9(tmpdir):
    testName = "rotational_group"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_10(tmpdir):
    testName = "switch_test"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh, deleteH=False)
    assert passGeo


def test_example_11(tmpdir):
    testName = "compact_bonding"
    thresh = 0.01
    passGeo = hp.runtestgeo_optonly(tmpdir, testName, thresh)
    assert passGeo


def test_example_12(tmpdir):
    testName = "triplebond_linear_broken"
    thresh = 0.01
    passGeo = hp.runtestgeo_optonly(tmpdir, testName, thresh)
    assert passGeo


def test_example_13(tmpdir):
    testName = "iodine_sulfur"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_14(tmpdir):
    testName = "oct_comp_greedy"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_15(tmpdir):
    testName = "atom_ordering_mismatch"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo


def test_example_16(tmpdir):
    testName = "iodide_radius"
    thresh = 0.01
    passGeo = hp.runtestgeo(tmpdir, testName, thresh)
    assert passGeo
