#!/usr/bin/env python
'''
    Copyright 2016 Kulik Lab @ MIT

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
# Dpt of Chemical Engineering, MIT

##########################################################
############  Main script that coordinates  ##############
#############  all parts of the program   ################
##########################################################

import sys, argparse, os, platform, shutil
from Scripts.inparse import *
from Scripts.generator import *

if __name__ == '__main__':
    ### run GUI by default ###
    args = sys.argv[1:]
    gui = True
    cmd = False
    try:
        import PyQt5
        from PyQt5.QtGui import *
        from Classes.mGUI import *
    except ImportError:
        if len(args)==0:
            print "\nGUI not supported since PyQt5 can not be loaded. Please use commandline version.\n"
            gui = False
    ####################################
    ### print help ###
    if '-h' in args or '-H' in args:
        # print help
        parser = argparse.ArgumentParser()
        args = parsecommandline(parser)
        parser.print_help()
        exit()
    ### run with gui ###
    elif gui and len(args)==0:
        ### create main application
        app = QApplication(sys.argv) # main application
        gui = mGUI(app) # main GUI class
        app.processEvents()
        app.exec_()
    ### if input file is specified run without GUI ###
    elif '-i' in args:
        gui = False
        # run from commandline
        emsg = startgen(sys.argv,False,gui)
    ### grab from commandline arguments ###
    else:
        gui = False
        infile = parseCLI(filter(None,args))
        args = ['main.py','-i',infile]
        emsg = startgen(args,False,gui)
        
        
        
        
        
        
        
        
        
