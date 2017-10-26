Instructions for Installing and Running Test Cases for molSimplify
      by Fang Liu  10/24/2017


I.  The structure of the repo

    |- molSimplify
    |      |- __init__.py
    |      |- Scripts
    |      |- Classes
    |      |_ ...
    |
    |- tests
    |      |-test_case1.py
    |      |-...
    |      |-test_caseN.py
    |       |- inputs
    |      |     |- case1.in
    |      |     |- ...
    |      |     |_ caseN.in
    |      |_ refs
    |            |- case1.xyz
    |            |- case1.report
    |            |_ ...
    |_ setup.py

    Reason for using this structure, please refer to: https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
    Notice that pytest automatically find test cases by look at python files
    named with
          test_XXX.py  XXX_test.py
    
    Under this stracture, the test files MUST have unique names

II. Add a test case
  1. Name the test input and output files in a consistent way:
     yourTestCase.in,  yourTestCase.xyz,  yourTestCase.report

  2. Put yourTestCase.in under: tests/inputs/
     Put yourTestCase.xyz,  yourTestCase.report  under: tests/refs/

  3. Create a test python script with the template shown below. Thresh
     is an optional tolerance for RMSD comparison, and defaults to 0.1 A.
     If you are adding a small test cas (e.g. hexachloride), consider
     reducing this parameter. Otherwise it is likely fine as is. 

############ test_yourTestCase.py  #########
import helperFuncs as hp

def test_example_1(tmpdir):
    [pass_xyz,pass_report] = hp.runtest(tmpdir,"yourTestCase",thresh)
    assert pass_xyz and pass_report

#############################################

III. Run Test Cases
  1. Tests will be automatically run by Travis. However, it's higly recomended
     that one runs the tests locally before commiting or push to repository to
     avoid breaking things

  2. Here are some different ways to run the tests. Assume that the we are at the
     root directory of the repository
     (1) make use of the setup.py:
         python setup.py test
     (2) Use the pytest command:
         pytest

     These will go through all test cases and report whether each is passed. If
     a test fails, pytest will automatically trace back to the function where
     it breaks. Standard output of molSimplify will not be printed unless the 
     test fails (then the output will be captured by pytest)

     If one wants to see all the stanford output on the screen, we can try:
     (1) Run pytest without capture:
         pytest -s

     (2) Manually run each test case under the tests/ directory
         python test_yourTestCase.py

