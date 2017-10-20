import pytest
import pybel, argparse
import os
import openbabel as ob
from molSimplify.Scripts.geometry import kabsch
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D
from pkg_resources import resource_filename, Requirement

def xyz2mol3D(xyz):
    mymol = mol3D()
    mymol.OBmol = pybel.readfile('xyz', xyz).next()
    mymol.convert2mol3D()
    mymol.OBmol = False
    return mymol

def fuzzy_compare_xyz(xyz1,xyz2,thresh):
    fuzzyEqual=False
    mol1 = xyz2mol3D(xyz1)
    mol2 = xyz2mol3D(xyz2)
    rmsd12 = mol1.rmsd(mol2)
    if rmsd12 < thresh:
        fuzzyEqual=True
    return fuzzyEqual

def jobname(infile):
    name=os.path.basename(infile)
    name=name.replace(".in","")
    return name

def jobdir(infile):
    name = jobname(infile)
    homedir = os.path.expanduser("~")
    mydir=homedir+'/Runs/'+name
    return mydir

def parse4test(infile,tmpdir):
    name = jobname(infile)
    f=tmpdir.join(os.path.basename(infile))
    newname = f.dirname+"/"+os.path.basename(infile)
    data=open(infile).readlines()
    newdata=""
    hasJobdir = False
    hasName = False
    for line in data:
        if not (("-jobdir" in line) or ("-name" in line)):
            newdata+=line
    newdata+="-jobdir "+name+"\n"
    newdata+="-name "+name+"\n"
    print newdata
    f.write(newdata)
    print "Input file parsed for test is located: ",newname
    return newname
