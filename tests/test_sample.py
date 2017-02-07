import pytest
import pybel
import openbabel as ob
from molSimplify.Classes.mol3D import *

xyz ="""3
Water
O         -0.63040        3.45320        1.15380
H          0.34240        3.45320        1.33890
H         -1.05330        3.45320        2.04920
"""

@pytest.fixture(scope="module")
def testmol():
    mymol = mol3D()
    mymol.OBmol = pybel.readstring("xyz",xyz)
    mymol.convert2mol3D()
    mymol.OBmol = False
    return mymol

def test_convert2OBmol(testmol):
    testmol.convert2OBmol()
    teststr = testmol.OBmol.write(format="mol",filename=None)
    answer = pybel.readstring("xyz",xyz).write(format="mol",filename=None)
    ### testing to see if all lines except for the title of the mol file are the same
    assert teststr.splitlines()[1:] == answer.splitlines()[1:]

