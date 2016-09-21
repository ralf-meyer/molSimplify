# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT
import os, sys
import glob, re, math, random, string, numpy, pybel
from math import pi
from geometry import *
from Classes.atom3D import *
from Classes.mol3D import*
from Classes.globalvars import globalvars
from operator import add

###############################
def write_periodic_mol3d_to_qe(mol,cell_vector,path):
        ## set global properties
        unique_atoms  = mol.getAtomTypes()
        globs = globalvars()
        print(globs)
        ## start writing this file:
        with open(path,'w') as f: 
                f.write("&CONTROL\n")
                f.write('calculation = "scf" \n')
                f.write('pseudo_dir = "./pseudo"\n')
                f.write('out_dir = "./scr"\n')
                f.write("/ \n")
 
        with open(path,'a') as f: 
                f.write("&SYSTEM\n")
                f.write("ibrav = 0 \n")
                f.write("celldm = \n")
                f.write("nat  = " + str(mol.natoms) + "\n")
                f.write("ntyp = " + str(len(unique_atoms)) + "\n")
                f.write('occupations = "smearing"\n')
                f.write('smearing = "gaussian"\n')
                f.write("/ \n")
        with open(path,'a') as f: 
                f.write("CELL_PARAMETERS {angstrom}\n")
                for cv in cell_vector:
                        ss = " ".join(str(e) for e in cv) + "\n"
                        f.write(ss) 
                f.write("/ \n")

        with open(path,'a') as f: 
                f.write("ATOMIC_SPECIES\n")
                for elements in unique_atoms:
                        if len(elements) == 1:
                               f.write(str(elements) + "     " + str(globs.amass()[elements][0]) + '\n') 
                        else:
                               f.write(str(elements) + "    " +  str(globs.amass()[elements][0]) + '\n') 
        with open(path,'a') as f: 
                f.write("ATOMIC_POSITIONS {angstrom}\n")
                for atom in mol.atoms:
                    xyz = atom.coords()
                    f.write("%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2]))
                f.write("/ \n")
                






