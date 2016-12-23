#!/usr/bin/env python
# Written by Terry for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
##########   Compiles predefined ligands in    ###########
##########   Ligands folder into a database    ###########
##########        for screening purposes       ###########
##########################################################

from Classes.globalvars import *
import argparse, sys, os, shutil, itertools, random
import pybel
globs = globalvars()

f = open(globs.installdir+'/Data/predefined_ligands.sdf','w')
ligfiles = os.listdir(globs.installdir+'/Ligands/')
for lig in ligfiles:
	if '.mol' in lig:
		g = open(globs.installdir+'/Ligands/'+lig,'r')
		s = g.read()
		g.close()
		ss = s.splitlines()
		for l in ss:
			f.write(l+'\n')
		f.write('\n$$$$\n')
f.close()		
