import pytest
import pybel, argparse
import os
import openbabel as ob
from molSimplify.Scripts.inparse import *
from molSimplify.Scripts.generator import *
from molSimplify.Classes.globalvars import *
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

def compare_report(report1,report2):
    data1=open(report1,'r').readlines()
    data2=open(report2,'r').readlines()
    Equal = True
    for i,lines in enumerate(data1):
        if Equal:
            Equal = (lines.strip() == data2[i].strip())
    return Equal

def runtest(tmpdir,name,thresh = 0.1):
    infile = resource_filename(Requirement.parse("molSimplify"),"tests/inputs/"+name+".in")
    newinfile = parse4test(infile,tmpdir)
    args =['main.py','-i', newinfile]
    startgen(args,False,False)
    myjobdir=jobdir(infile)
    output_xyz = myjobdir + '/'+ name + '.xyz'
    output_report = myjobdir + '/'+ name + '.report'
    ref_xyz = resource_filename(Requirement.parse("molSimplify"),"tests/refs/"+name+".xyz")
    ref_report = resource_filename(Requirement.parse("molSimplify"),"tests/refs/"+name+".report")
    print "Test input file: ", newinfile
    print "Test output files are generated in ",myjobdir
    print "Output xyz file: ", output_xyz
    pass_xyz=fuzzy_compare_xyz(output_xyz,ref_xyz,thresh)
    pass_report = compare_report(output_report,ref_report)
    print "Reference xyz file: ", ref_xyz
    print "Test report file: ", output_report
    print "Reference report file: ", ref_report
    print "Reference xyz status: ", pass_xyz
    print "Reference report status: ", pass_report
    return [pass_xyz,pass_report]
