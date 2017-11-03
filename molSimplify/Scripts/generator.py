#!/usr/bin/env python
'''
    Copyright 2017 Kulik Lab @ MIT

    This file is part of molSimplify.
    molSimplify is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    molSimplify is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with molSimplify. If not, see http://www.gnu.org/licenses/.
'''

# Written by Tim Ioannidis for HJK Group
# Extended by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
############  Main script that coordinates  ##############
#############  all parts of the program   ################
##########################################################

import sys, os, random, shutil, inspect, argparse, openbabel
from molSimplify.Scripts.rungen import *
from molSimplify.Scripts.io import *
from molSimplify.Scripts.inparse import *
from molSimplify.Scripts.dbinteract import *
from molSimplify.Scripts.postproc import *
from molSimplify.Scripts.cellbuilder import*
from molSimplify.Scripts.chains import*
from molSimplify.Scripts.findcorrelations import*
from molSimplify.Classes.globalvars import *
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D
from math import sqrt
from math import floor


def startgen(argv,flag,gui):
    emsg = False
    ### check for configuration file ##
    homedir = os.path.expanduser("~")
    #configfile = False if not glob.glob(homedir+'/.molSimplify') else True
    #if not configfile:
    #    print "It looks like the configuration file '~/.molSimplify' does not exist!Please follow the next steps to configure the file."
    #    instdir = raw_input("Please select the full path of the top installation directory for the program: ")
    #    cdbdir = raw_input("Please specify the full path of the directory containing chemical databases:")
    #    mwfn = raw_input("Specify the full path to the Multiwfn executable (for post-processing):")
    #    f = open(homedir+'/.molSimplify','w')
    #    if len(instdir) > 1:
    #        f.write("INSTALLDIR="+instdir+'\n')
    #    if len(cdbdir) > 1:
    #        f.write("CHEMDBDIR="+cdbdir+'\n')
    #    if len(mwfn) > 1 :
    #        f.write("MULTIWFN="+mwfn[0]+'\n')
    #    f.close()
    ### end set-up configuration file ###
    ############ GLOBALS DEFINITION ############
    globs = globalvars()
    #installdir = globs.installdir
    rundir = globs.rundir
    PROGRAM = globs.PROGRAM
    ###### END GLOBALS DEFINITION ##############
    # correct installdir
    #if installdir[-1]!='/':
    #    installdir+='/'
    # print welcome message
    ss = "\n************************************************************"
    ss += "\n******** Welcome to "+PROGRAM+"! Let's get started. ********\n"
    ss += "************************************************************\n\n"
    if not flag:
        print ss
    sys.argv = argv
    parser = argparse.ArgumentParser()
    args = parseall(parser)
    # check if input file exists
    if not glob.glob(args.i):
        emsg = 'Input file '+args.i+' does not exist. Please specify a valid input file.\n'
        print emsg
        return emsg
    args.gui = gui # add gui flag
        # parse input file
    if args.i:
        parseinputfile(args)
    if not args.postp and not args.dbsearch and not args.dbfinger and not args.drawmode and not (args.slab_gen or args.place_on_slab) and not (args.chain) and not (args.correlate): # check input arguments
        # check input arguments
        print 'Checking input...'
        if args.tsgen:
            emsg = checkinput_ts(args)
        else:
            emsg = checkinput(args)
        # check before cleaning input arguments and clean only if checked
        cleaninput(args)
    args.gui = False # deepcopy will give error
    if emsg:
        del args
        return emsg
    # check for jobs directory
    rundir = args.rundir+'/' if (args.rundir) else rundir
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    ################### START MAIN ####################
    args0 = copy.deepcopy(args) # save initial arguments
    # add gui flag
    args.gui = gui
    # postprocessing run?

    if (args.postp):
        postproc(rundir,args,globs)
    # database search?
    elif (args.dbsearch or args.dbfinger):
        emsg = dbsearch(rundir,args,globs)
        if emsg:
            del args
            return emsg
        else:
            print 'Successful database search!\n'
    # random generation?
    elif (args.rgen): # check if random generation was requested
        if args.charge:
            args.charge = args.charge[0]
        if args.spin:
            args.spin = args.spin[0]
        corests=args.core
        for cc in corests:
            args = copy.deepcopy(args0)
            # add gui flag
            args.gui = gui
            args.core = cc
            if (args.lig or args.coord or args.lignum or args.ligocc): # constraints given?
                args, emsg = constrgen(rundir,args,globs)
                if emsg:
                    del args
                    return emsg
            else:
                emsg = 'For random generation specify at least a ligand, coordination or ligand types.\n'
                print emsg
                del args
                return emsg
    elif args.drawmode:
        emsg = draw_supervisor(args,rundir)            
    # slab/place on slab?
    elif (args.slab_gen or args.place_on_slab):
        emsg = slab_module_supervisor(args,rundir)
    # chain builder
    elif (args.chain):
        print('chain on')
        emsg = chain_builder_supervisor(args,rundir)
    # correlation analysis
    elif (args.correlate):
        print('analysis is looking for correlations')
        analysis_supervisor(args,rundir)
    # normal structure generation or transition state building
    else:
        args = copy.deepcopy(args0)
        # add gui flag
        args.gui = gui
        corests=args.core
        if args.tsgen: # goes through multigenruns for maximum interoperability
            print('building a transition state')
        else:
            print('building an equilibrium complex')
        for cc in corests:
            args.core = cc
            emsg = multigenruns(rundir,args,globs)
            if emsg:
                print emsg
                del args
                return emsg
    ss =  "\n**************************************************************"
    ss += "\n***** Thank you for using "+PROGRAM+". Have a nice day! ******\n"
    ss += "**************************************************************"
    ss += globs.about
    if not flag:
        print ss
    del args
    return emsg
    #################### END MAIN #####################

if __name__ == "__main__":
    startgen()
