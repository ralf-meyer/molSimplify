import pytest
import pybel, argparse
import os
import openbabel as ob
import helperFuncs as hp
from molSimplify.Scripts.inparse import *
from molSimplify.Scripts.generator import *
from molSimplify.Classes.globalvars import *
from pkg_resources import resource_filename, Requirement

def test_example1(tmpdir):
    infile = resource_filename(Requirement.parse("molSimplify"),"tests/inputs/example-1.in")
    name = hp.jobname(infile)
    newinfile = hp.parse4test(infile,tmpdir)
    args =['main.py','-i', newinfile]
    startgen(args,False,False)
    myjobdir=hp.jobdir(infile)
    output_xyz = myjobdir + '/'+ name + '.xyz'
    ref_xyz = resource_filename(Requirement.parse("molSimplify"),"tests/refs/"+name+".xyz")
    print "Test output files are generated in ",myjobdir
    print "Output xyz file: ", output_xyz
    print "Reference xyz file: ", ref_xyz
    assert hp.fuzzy_compare_xyz(output_xyz,ref_xyz,0.001)
