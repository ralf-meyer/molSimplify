# Written by JP for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
##########  Top level script that coordinates  ###########
##########    generation of structues, input   ###########
##########         files, jobscripts           ###########
##########################################################

from structgen import *
from io import *
import argparse, sys, os, shutil, itertools, random
import pybel

def name_complex(core,ligs,ligoc,args):
    center = core.getAtom(0).symbol()
    name = center + '_'
    if args.oxstate:
        ox = str(args.oxstate)
    else:
        ox = "0"
    name += "_ " + str(ox)
    if args.spin:
        spin = str(args.spin)
    else:
        spin = "0"
    name += "_ " + str(spin)
    for i,lig in enumerate(ligs):
        names += '_' + str(lig[:3]) + '-' + str(ligoc[i])
    names += "_"+str(spin)
    return name

