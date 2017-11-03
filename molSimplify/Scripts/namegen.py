# Written by JP for HJK Group
# Dpt of Chemical Engineering, MIT

from structgen import *
from molSimplify.Scripts.io import *
import argparse, sys, os, shutil, itertools, random

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
