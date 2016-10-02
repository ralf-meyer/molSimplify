# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

####################################################
########### Defines main class of GUI  #############
####################################################
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from Classes.mWidgets import *
from Classes.globalvars import *
from Classes.mol3D import mol3D
from Scripts.generator import startgen
from Scripts.grabguivars import *
from Scripts.io import *
from Scripts.addtodb import *
import sys, os, random, shutil, unicodedata, inspect, glob, time
import pybel

class mGUI():
    getgeoms()
    ### constructor of gui ###
    def __init__(self,app):
        # build gui
        self.app = app
        self.initGUI(app)

    ### builds the gui
    def initGUI(self,app):
        '''
        ######################
        ### build main GUI ###
        ######################
        '''
        ### check for configuration file ###
        homedir = os.path.expanduser("~")
        globs = globalvars() # global variables
        overX = True if 'localhost' in os.environ['DISPLAY'].lower() else False # detect running over X
        configfile = False if not glob.glob(homedir+'/.molSimplify') else True
        if not configfile:
            self.wwindow = mQMainWindow() 
            self.wwindow.resize(0.5,0.5)
            QMessageBox.information(self.wwindow,'Setup',"It looks like the configuration file '~/.molSimplify' does not exist!Please follow the next steps to configure the file.")
            QMessageBox.information(self.wwindow,'Installation directory',"Please select the top installation directory for the program.")
            instdir = QFileDialog.getExistingDirectory(self.wwindow, "Select Directory")
            f = open(homedir+'/.molSimplify','w')
            if len(instdir) > 1: 
                f.write("INSTALLDIR="+instdir+'\n')
        ### end set-up configuration file ###
        ### main window widget
        self.wmwindow = mQMainWindow()
        self.wmwindow.show()
        self.wmain = QWidget()
        self.wmwindow.setWindowTitle("molSimplify")
        self.wmwindow.setMinimumSize(1000,700)
        # set background color
        p = QPalette()
        p.setColor(QPalette.Background,Qt.white)
        self.wmwindow.setPalette(p)
        ### main grid layout ###
        self.grid = QGridLayout()
        ### stacked layouts ###
        self.sgrid = QStackedLayout()
        self.sgrid.setStackingMode(1)
        self.sgrid.addWidget(self.wmain)
        ### create menubar and callbacks ###
        menubar = self.wmwindow.menuBar()
        menu0 = menubar.addMenu('&File')
        menu1 = menubar.addMenu('&Load')
        menu2 = menubar.addMenu('&Help')
        exitAction = QAction('&Exit',self.wmwindow) 
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.qexit)
        menu0.addAction(exitAction)
        saveAction = QAction('&Save As..',self.wmwindow)
        saveAction.setShortcut('Ctrl+S')
        saveAction.setStatusTip('Save current input settings')
        saveAction.triggered.connect(self.qsaveinput)
        menu0.addAction(saveAction)
        loadAction = QAction('&Load',self.wmwindow)
        loadAction.setShortcut('Ctrl+O')
        loadAction.setStatusTip('Load input file')
        loadAction.triggered.connect(self.qloadinput)
        menu1.addAction(loadAction)
        helpAction = QAction('&About',self.wmwindow)
        helpAction.setShortcut('Ctrl+H')
        helpAction.setStatusTip('Show program information')
        helpAction.triggered.connect(self.qshowhelp)
        menu2.addAction(helpAction)
        ### place title top ###
        self.grid.setRowMinimumHeight(0,15)
        self.grid.setRowMinimumHeight(2,120)
        self.grid.setRowMinimumHeight(3,15)
        self.grid.setRowMinimumHeight(4,50)
        clogo = mQPixmap(globs.installdir+'/icons/logo.png')
        self.grid.addWidget(clogo,1,9,2,18)
        self.txtdev = mQLabel('Developed by Kulik group @ MIT','','c',16)
        self.grid.addWidget(self.txtdev,19,9,2,18)
        ###################################################
        ###################################################
        ########## STRUCTURE GENERATION INPUTS ############
        ###################################################
        ###################################################
        ### title ###
        self.txtgpar = mQLabel('Structure specification','','c',20)
        self.grid.addWidget(self.txtgpar,4,0,2,14)
        ### core structure specification ###
        ctip = 'Core of structure'
        self.rtcore = mQLabel('Core:',ctip,'C',12)
        f = QFont("Helvetica",14,75)
        self.rtcore.setFont(f)
        self.grid.addWidget(self.rtcore,7,2,1,1)
        self.etcore = mQLineEdit('',ctip,'l',12)
        self.grid.addWidget(self.etcore,7,3,1,3)
        ### Connection atoms for core ###
        ctip = 'Specify connection atoms for core if using SMILES or custom cores (default: 1)'
        self.rtccat = mQLabel('Core connections:',ctip,'C',12)
        self.rtccat.setWordWrap(True)
        self.grid.addWidget(self.rtccat,7,6,1,2)
        self.etccat = mQLineEdit('',ctip,'l',12)
        self.grid.addWidget(self.etccat,7,8,1,2)
        ### replace option ###
        ctip = 'Replace ligand at specified connection point'
        self.replig = mQCheckBox('replace',ctip,12)
        self.grid.addWidget(self.replig,7,10,1,4)
        # coordination
        ctip = 'Number of ligands connected to the metal.'
        self.rcoord = mQLabel('Coordination:',ctip,'Cr',12)
        self.grid.addWidget(self.rcoord,8,2,1,2)
        coords,geomnames,geomshorts,geomgroups = getgeoms()
        qcav = sorted(list(set(coords)))
        self.dcoord = mQComboBox(qcav,ctip,12)
        self.dcoord.setCurrentIndex(5)
        self.dcoord.currentIndexChanged.connect(self.matchgeomcoord)
        self.grid.addWidget(self.dcoord,8,4,1,1)
        # geometry of coordination
        self.dcoordg = mQComboBox('','',12)
        self.dcoordg.setCurrentIndex(0)
        self.matchgeomcoord()
        self.grid.addWidget(self.dcoordg,8,5,1,3)
        # add new coordination
        ctip = 'Add geometry'
        self.butaddg = mQPushButton('Add geometry',ctip,12)
        self.butaddg.clicked.connect(self.addgeom)
        self.grid.addWidget(self.butaddg,8,8,1,2)
        # add new coordination
        ctip = 'View geometry with atom labels'
        self.butvg = mQPushButton('View geometry',ctip,12)
        self.butvg.clicked.connect(self.viewgeom)
        self.grid.addWidget(self.butvg,8,10,1,1)
        #############################################
        ################## LIGANDS ##################
        ### ligands tables ###
        ctip0 = 'Ligand(s) to be used' 
        self.rtligh = mQLabel('Ligand',ctip0,'c',12) # ligand header
        f = QFont("Helvetica",12,75)
        self.rtligh.setFont(f)
        ctip1 = 'Occurrence of corresponding ligand(s)'
        self.rtligocch = mQLabel('Frequency',ctip1,'c',12) # occurrence header
        ctip2 = 'Connection atom(s) of ligands (default: 1).'
        self.rtsmicath = mQLabel('Connections',ctip2,'c',12) # connection atom header
        ctip3 = 'Do not remove hydrogens while connecting ligand to core. default False' # keep Hs header
        self.keepHh = mQLabel('keep\nHs',ctip3,'c',12) # occurrence header
        ctip4 = 'Custom bond length for M-L in Angstrom' 
        self.MLbondsh = mQLabel('M-L\nbond',ctip4,'c',12)# custom metal ligand bond length header
        ctip5 = 'Custom angles for connection points (polar theta, azimuthal phi) in degrees separated with /. e.g 10/20'
        self.canglesh = mQLabel('Angle',ctip5,'c',12) # custom angles for distortion header
        ctip6 = 'Name of ligand'
        self.nameligh = mQLabel('Name',ctip6,'c',12) # name of ligand header
        ctip7 = 'Force ligand order and disable smart reordering'
        self.ligfloc = mQCheckBox('Force\nlocation',ctip7,12) 
        ctip8 = 'Ligand smart alignment. Aligns first the bulky ligands.'
        self.ligfalign = mQCheckBox('Smart\nalignment',ctip8,12) 
        self.ligfalign.setChecked(True)
        # add to layout
        self.grid.setColumnMinimumWidth(0,150)
        self.grid.addWidget(self.rtligh,9,0,1,2)
        self.grid.addWidget(self.rtligocch,9,2,1,2)
        self.grid.addWidget(self.rtsmicath,9,4,1,2)
        self.grid.addWidget(self.keepHh,9,6,1,1)
        self.grid.addWidget(self.MLbondsh,9,7,1,2)
        self.grid.addWidget(self.canglesh,9,9,1,1)
        self.grid.addWidget(self.nameligh,9,10,1,1)
        self.grid.addWidget(self.ligfalign,8,13,1,1)
        self.grid.addWidget(self.ligfloc,9,13,1,1)
        ## ligands ##
        self.lig0 = mQLineEdit('',ctip0,'l',12)
        self.lig1 = mQLineEdit('',ctip0,'l',12)
        self.lig2 = mQLineEdit('',ctip0,'l',12)
        self.lig3 = mQLineEdit('',ctip0,'l',12)
        self.lig4 = mQLineEdit('',ctip0,'l',12)
        self.lig5 = mQLineEdit('',ctip0,'l',12)
        self.lig6 = mQLineEdit('',ctip0,'l',12)
        self.lig7 = mQLineEdit('',ctip0,'l',12)
        # disable
        self.lig1.setDisabled(True)
        self.lig2.setDisabled(True)
        self.lig3.setDisabled(True)
        self.lig4.setDisabled(True)
        self.lig5.setDisabled(True)
        self.lig6.setDisabled(True)
        self.lig7.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0,10,0,1,2)
        self.grid.addWidget(self.lig1,11,0,1,2)
        self.grid.addWidget(self.lig2,12,0,1,2)
        self.grid.addWidget(self.lig3,13,0,1,2)
        self.grid.addWidget(self.lig4,14,0,1,2)
        self.grid.addWidget(self.lig5,15,0,1,2)
        self.grid.addWidget(self.lig6,16,0,1,2)
        self.grid.addWidget(self.lig7,17,0,1,2)
        ## occurrences ##
        self.lig0occ = mQSpinBox(ctip1)
        self.lig1occ = mQSpinBox(ctip1)
        self.lig2occ = mQSpinBox(ctip1)
        self.lig3occ = mQSpinBox(ctip1)
        self.lig4occ = mQSpinBox(ctip1)
        self.lig5occ = mQSpinBox(ctip1)
        self.lig6occ = mQSpinBox(ctip1)
        self.lig7occ = mQSpinBox(ctip1)
        # disable
        self.lig1occ.setDisabled(True)
        self.lig2occ.setDisabled(True)
        self.lig3occ.setDisabled(True)
        self.lig4occ.setDisabled(True)
        self.lig5occ.setDisabled(True)
        self.lig6occ.setDisabled(True)
        self.lig7occ.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0occ,10,2,1,2)
        self.grid.addWidget(self.lig1occ,11,2,1,2)
        self.grid.addWidget(self.lig2occ,12,2,1,2)
        self.grid.addWidget(self.lig3occ,13,2,1,2)
        self.grid.addWidget(self.lig4occ,14,2,1,2)
        self.grid.addWidget(self.lig5occ,15,2,1,2)
        self.grid.addWidget(self.lig6occ,16,2,1,2)
        self.grid.addWidget(self.lig7occ,17,2,1,2)
        ## connections ##
        self.lig0conn = mQLineEdit('',ctip2,'l',12)
        self.lig1conn = mQLineEdit('',ctip2,'l',12)
        self.lig2conn = mQLineEdit('',ctip2,'l',12)
        self.lig3conn = mQLineEdit('',ctip2,'l',12)
        self.lig4conn = mQLineEdit('',ctip2,'l',12)
        self.lig5conn = mQLineEdit('',ctip2,'l',12)
        self.lig6conn = mQLineEdit('',ctip2,'l',12)
        self.lig7conn = mQLineEdit('',ctip2,'l',12)
        # disable
        self.lig1conn.setDisabled(True)
        self.lig2conn.setDisabled(True)
        self.lig3conn.setDisabled(True)
        self.lig4conn.setDisabled(True)
        self.lig5conn.setDisabled(True)
        self.lig6conn.setDisabled(True)
        self.lig7conn.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0conn,10,4,1,2)
        self.grid.addWidget(self.lig1conn,11,4,1,2)
        self.grid.addWidget(self.lig2conn,12,4,1,2)
        self.grid.addWidget(self.lig3conn,13,4,1,2)
        self.grid.addWidget(self.lig4conn,14,4,1,2)
        self.grid.addWidget(self.lig5conn,15,4,1,2)
        self.grid.addWidget(self.lig6conn,16,4,1,2)
        self.grid.addWidget(self.lig7conn,17,4,1,2)
        ## keep Hydrogens #
        self.lig0H = mQComboBox(['no','yes'],ctip3,12)
        self.lig1H = mQComboBox(['no','yes'],ctip3,12)
        self.lig2H = mQComboBox(['no','yes'],ctip3,12)
        self.lig3H = mQComboBox(['no','yes'],ctip3,12)
        self.lig4H = mQComboBox(['no','yes'],ctip3,12)
        self.lig5H = mQComboBox(['no','yes'],ctip3,12)
        self.lig6H = mQComboBox(['no','yes'],ctip3,12)
        self.lig7H = mQComboBox(['no','yes'],ctip3,12)
        # disable
        self.lig1H.setDisabled(True)
        self.lig2H.setDisabled(True)
        self.lig3H.setDisabled(True)
        self.lig4H.setDisabled(True)
        self.lig5H.setDisabled(True)
        self.lig6H.setDisabled(True)
        self.lig7H.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0H,10,6,1,1)
        self.grid.addWidget(self.lig1H,11,6,1,1)
        self.grid.addWidget(self.lig2H,12,6,1,1)
        self.grid.addWidget(self.lig3H,13,6,1,1)
        self.grid.addWidget(self.lig4H,14,6,1,1)
        self.grid.addWidget(self.lig5H,15,6,1,1)
        self.grid.addWidget(self.lig6H,16,6,1,1)
        self.grid.addWidget(self.lig7H,17,6,1,1)
        ## ML bonds ##
        self.lig0ML = mQLineEdit('',ctip4,'l',12)
        self.lig1ML = mQLineEdit('',ctip4,'l',12)
        self.lig2ML = mQLineEdit('',ctip4,'l',12)
        self.lig3ML = mQLineEdit('',ctip4,'l',12)
        self.lig4ML = mQLineEdit('',ctip4,'l',12)
        self.lig5ML = mQLineEdit('',ctip4,'l',12)
        self.lig6ML = mQLineEdit('',ctip4,'l',12)
        self.lig7ML = mQLineEdit('',ctip4,'l',12)
        # disable
        self.lig1ML.setDisabled(True)
        self.lig2ML.setDisabled(True)
        self.lig3ML.setDisabled(True)
        self.lig4ML.setDisabled(True)
        self.lig5ML.setDisabled(True)
        self.lig6ML.setDisabled(True)
        self.lig7ML.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0ML,10,7,1,2)
        self.grid.addWidget(self.lig1ML,11,7,1,2)
        self.grid.addWidget(self.lig2ML,12,7,1,2)
        self.grid.addWidget(self.lig3ML,13,7,1,2)
        self.grid.addWidget(self.lig4ML,14,7,1,2)
        self.grid.addWidget(self.lig5ML,15,7,1,2)
        self.grid.addWidget(self.lig6ML,16,7,1,2)
        self.grid.addWidget(self.lig7ML,17,7,1,2)
        ## custom angles ##
        self.lig0an = mQLineEdit('',ctip5,'l',12)
        self.lig1an = mQLineEdit('',ctip5,'l',12)
        self.lig2an = mQLineEdit('',ctip5,'l',12)
        self.lig3an = mQLineEdit('',ctip5,'l',12)
        self.lig4an = mQLineEdit('',ctip5,'l',12)
        self.lig5an = mQLineEdit('',ctip5,'l',12)
        self.lig6an = mQLineEdit('',ctip5,'l',12)
        self.lig7an = mQLineEdit('',ctip5,'l',12)
        # disable
        self.lig1an.setDisabled(True)
        self.lig2an.setDisabled(True)
        self.lig3an.setDisabled(True)
        self.lig4an.setDisabled(True)
        self.lig5an.setDisabled(True)
        self.lig6an.setDisabled(True)
        self.lig7an.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0an,10,9,1,1)
        self.grid.addWidget(self.lig1an,11,9,1,1)
        self.grid.addWidget(self.lig2an,12,9,1,1)
        self.grid.addWidget(self.lig3an,13,9,1,1)
        self.grid.addWidget(self.lig4an,14,9,1,1)
        self.grid.addWidget(self.lig5an,15,9,1,1)
        self.grid.addWidget(self.lig6an,16,9,1,1)
        self.grid.addWidget(self.lig7an,17,9,1,1)
        ## ligand names ##
        self.lig0nam = mQLineEdit('',ctip6,'l',12)
        self.lig1nam = mQLineEdit('',ctip6,'l',12)
        self.lig2nam = mQLineEdit('',ctip6,'l',12)
        self.lig3nam = mQLineEdit('',ctip6,'l',12)
        self.lig4nam = mQLineEdit('',ctip6,'l',12)
        self.lig5nam = mQLineEdit('',ctip6,'l',12)
        self.lig6nam = mQLineEdit('',ctip6,'l',12)
        self.lig7nam = mQLineEdit('',ctip6,'l',12)
        # disable
        self.lig1nam.setDisabled(True)
        self.lig2nam.setDisabled(True)
        self.lig3nam.setDisabled(True)
        self.lig4nam.setDisabled(True)
        self.lig5nam.setDisabled(True)
        self.lig6nam.setDisabled(True)
        self.lig7nam.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0nam,10,10,1,1)
        self.grid.addWidget(self.lig1nam,11,10,1,1)
        self.grid.addWidget(self.lig2nam,12,10,1,1)
        self.grid.addWidget(self.lig3nam,13,10,1,1)
        self.grid.addWidget(self.lig4nam,14,10,1,1)
        self.grid.addWidget(self.lig5nam,15,10,1,1)
        self.grid.addWidget(self.lig6nam,16,10,1,1)
        self.grid.addWidget(self.lig7nam,17,10,1,1)
        ## add buttons ##
        ctip = 'Add new ligand.'
        self.lig0add = mQPushButton('+',ctip,12)
        self.lig1add = mQPushButton('+',ctip,12)
        self.lig2add = mQPushButton('+',ctip,12)
        self.lig3add = mQPushButton('+',ctip,12)
        self.lig4add = mQPushButton('+',ctip,12)
        self.lig5add = mQPushButton('+',ctip,12)
        self.lig6add = mQPushButton('+',ctip,12)
        # add callbacks
        self.lig0add.clicked.connect(self.addlig0)
        self.lig1add.clicked.connect(self.addlig1)
        self.lig2add.clicked.connect(self.addlig2)
        self.lig3add.clicked.connect(self.addlig3)
        self.lig4add.clicked.connect(self.addlig4)
        self.lig5add.clicked.connect(self.addlig5)
        self.lig6add.clicked.connect(self.addlig6)
        # hide and disable
        self.lig1add.hide()
        self.lig1add.setDisabled(True)
        self.lig2add.hide()
        self.lig2add.setDisabled(True)
        self.lig3add.hide()
        self.lig3add.setDisabled(True)
        self.lig4add.hide()
        self.lig4add.setDisabled(True)
        self.lig5add.hide()
        self.lig5add.setDisabled(True)
        self.lig6add.hide()
        self.lig6add.setDisabled(True)
        # add to layout
        self.grid.addWidget(self.lig0add,10,13,1,1)
        self.grid.addWidget(self.lig1add,11,13,1,1)
        self.grid.addWidget(self.lig2add,12,13,1,1)
        self.grid.addWidget(self.lig3add,13,13,1,1)
        self.grid.addWidget(self.lig4add,14,13,1,1)
        self.grid.addWidget(self.lig5add,15,13,1,1)
        self.grid.addWidget(self.lig6add,16,13,1,1)
        #############################################
        ## Draw ligand button ##
        ctip = 'Generate 2D ligand representation.'
        self.butDrl = mQPushButton('Draw ligands',ctip,12)
        self.butDrl.clicked.connect(self.drawligs)
        self.grid.addWidget(self.butDrl,18,2,1,2)
        ## Search DB button ##
        ctip = 'Search for ligands in chemical databases.'
        self.searchDB = mQPushButton('Search DB',ctip,12)
        self.searchDB.clicked.connect(self.searchDBW)
        self.grid.addWidget(self.searchDB,18,4,1,3)
        ## Local database button ##
        ctip = 'Add core/ligand/binding species to local database.'
        self.butADB = mQPushButton('Add to local DB',ctip,12)
        self.butADB.clicked.connect(self.enableDB)
        self.grid.addWidget(self.butADB,18,7,1,3)
        ##################################################
        ##################################################
        ########### GENERAL PARAMETERS INPUTS ############
        ##################################################
        ##################################################
        self.txtgp = mQLabel('General parameters','','c',20)
        self.grid.addWidget(self.txtgp,4,21,2,9)
        # random generation
        ctip = 'Enable random generation.'
        self.randomchk = mQCheckBox('Random generation',ctip,12)
        self.randomchk.stateChanged.connect(self.enablerandom)
        self.grid.addWidget(self.randomchk,7,23,1,2)
        # charge calculation
        ctip = 'Calculate charge based on ox state and ligands'
        self.chch = mQCheckBox('Calculate charge',ctip,12)
        self.chch.setDisabled(True)
        self.grid.addWidget(self.chch,7,25,1,2)
        # number of random generated structures
        ctip = 'Number of structures to be randomly generated.'
        self.rtrgen = mQLabel('Structures:',ctip,'Cr',12)
        self.etrgen = mQLineEdit('',ctip,'l',12)
        self.grid.addWidget(self.rtrgen,8,23,1,1)
        self.grid.addWidget(self.etrgen,8,24,1,1)
        self.rtrgen.setDisabled(True)
        self.etrgen.setDisabled(True)
        # number of different ligands to use
        ctip = 'For random generation: total number of different ligands including specified ones.'
        self.rtlignum = mQLabel('Different\nligands:',ctip,'Cr',12)
        qcav = ['1','2','3','4','5','6','7','8','9','10']
        self.etlignum = mQComboBox(qcav,ctip,12)
        self.grid.addWidget(self.rtlignum,8,25,1,1)
        self.grid.addWidget(self.etlignum,8,26,1,1)
        self.rtlignum.setDisabled(True)
        self.etlignum.setDisabled(True)
        # keep Hs
        ctip = 'Keep hydrogens for random generation.'
        self.randkHs = mQCheckBox('Keep Hs',ctip,12)
        self.grid.addWidget(self.randkHs,8,27,1,1)
        self.randkHs.setDisabled(True)
        # group of different ligands
        ctip = 'For random generation: select random ligands from group.'
        self.rtliggrp = mQLabel('Ligand\ngroup:',ctip,'Cr',12)
        qcav0 = getligroups(readdict(globs.installdir+'/Ligands/ligands.dict'))
        qcav = filter(None,qcav0.split(' '))
        self.etliggrp = mQComboBox(qcav,ctip,12)
        self.grid.addWidget(self.rtliggrp,9,23,1,1)
        self.grid.addWidget(self.etliggrp,9,24,1,1)
        self.rtliggrp.setDisabled(True)
        self.etliggrp.setDisabled(True)
        # ligand category
        ctip = 'For random generation: select random ligands from category.\nOptions are all'
        ctip += '"build" for building complexes, "functinoalize" for functionalizing'
        self.rtligctg = mQLabel('Ligand\ncategory:',ctip,'Cr',12)
        qcav = ['all','build','functionalize']
        self.etligctg = mQComboBox(qcav,ctip,12)
        self.grid.addWidget(self.rtligctg,9,25,1,1)
        self.grid.addWidget(self.etligctg,9,26,1,2)
        self.rtligctg.setDisabled(True)
        self.etligctg.setDisabled(True)
        # oxidation state
        ctip = 'Metal Oxidation state'
        self.roxstate = mQLabel('Ox State:',ctip,'Cr',12)
        qcav = ['0','I','II','III','IV','V','VI','VII','VIII']
        self.doxs = mQComboBox(qcav,ctip,12)
        self.doxs.setCurrentIndex(0)
        self.grid.addWidget(self.roxstate,10,23,1,1)
        self.grid.addWidget(self.doxs,10,24,1,1)
        # spin state
        ctip = 'System spin multiplicity'
        self.rspstate = mQLabel('Spin:',ctip,'Cr',12)
        qcav = ['1','2','3','4','5','6','7','8','9','10']
        self.dspin = mQComboBox(qcav,ctip,12)
        self.dspin.setCurrentIndex(0)
        self.grid.addWidget(self.rspstate,10,26,1,1)
        self.grid.addWidget(self.dspin,10,27,1,1)
        # force field optimization
        ctip = 'Perform Force Field optimization'
        self.chkFF = mQCheckBox('FF optimize',ctip,12)
        self.chkFF.stateChanged.connect(self.enableffinput)
        self.grid.addWidget(self.chkFF,11,23,1,2)
        # generate all
        ctip = 'Generate structure with and without optimization.'
        self.chkgenall = mQCheckBox('Generate all',ctip,12)
        self.chkgenall.stateChanged.connect(self.disableffinput)
        self.grid.addWidget(self.chkgenall,12,23,1,2)
        # perform optimization
        ctip = 'Select Force Field'
        qcav = ['MMFF94','UFF','gchemical','GAFF']
        self.dff = mQComboBox(qcav,ctip,12)
        self.dff.setCurrentIndex(0)
        self.dff.setDisabled(True)
        self.grid.addWidget(self.dff,11,25,1,3)
        # optimize before or after
        ctip = 'Optimize before or after building the structure'
        qcav = ['Before','After','Before & After']
        self.dffba = mQComboBox(qcav,ctip,12)
        self.dffba.setDisabled(True)
        self.dffba.setCurrentIndex(2)
        self.grid.addWidget(self.dffba,12,25,1,3)
        # create distortion slider
        ctip = 'Percent random distortion from default coordination geometry.'
        self.distper = mQLabel('Distort:0%',ctip,'Cr',12)
        self.sdist = mQSlider(ctip)
        self.sdist.valueChanged.connect(self.sliderChanged)
        self.grid.addWidget(self.distper,13,24,1,2)
        self.grid.addWidget(self.sdist,13,26,1,2)
        ### jobs dir ###
        ctip = 'Top directory for job folders.'
        self.rtrdir = mQLabel('Jobs dir:',ctip,'Cr',12)
        self.grid.addWidget(self.rtrdir,14,23,1,1)
        self.etrdir = mQLineEdit(globs.rundir,ctip,'l',12)
        self.grid.addWidget(self.etrdir,14,24,1,2)
        # button for browsing rundir
        ctip = 'Browse running directory.'
        self.butpbrdir = mQPushButton('Browse..',ctip,12)
        self.butpbrdir.clicked.connect(self.dirload)
        self.grid.addWidget(self.butpbrdir,14,26,1,2)
        # suffix
        ctip = 'Suffix for job directories.'
        self.rtsuff = mQLabel('Suffix:',ctip,'Cr',12)
        self.etsuff = mQLineEdit('',ctip,'l',12)
        self.grid.addWidget(self.rtsuff,15,23,1,1)
        self.grid.addWidget(self.etsuff,15,24,1,2)
        # structure generation
        ctip = 'Generate structures'
        self.butGen = mQPushButton('Generate',ctip,18)
        self.butGen.clicked.connect(self.runGUI)
        self.grid.addWidget(self.butGen,17,23,2,2)
        # post-processing setup
        ctip = 'Setup post-processing'
        self.butPost = mQPushButton('Post-process',ctip,16)
        self.butPost.clicked.connect(self.setupp)
        self.grid.addWidget(self.butPost,17,26,2,2)
        ###################################################
        ###################################################
        ########### ADDITIONAL MOLECULE INPUTS ############
        ###################################################
        ###################################################
        ### generate edit texts for additional molecule
        self.txtamol = mQLabel('Additional molecule','','c',20)
        self.txtamol.setDisabled(True)
        self.grid.addWidget(self.txtamol,4,34,2,5)
        # additional molecule
        ctip = 'Place additional molecule'
        self.chkM = mQCheckBox('Extra molecule',ctip,12)
        self.chkM.stateChanged.connect(self.enableemol)
        self.grid.addWidget(self.chkM,7,35,1,3)
        # name of binding species
        ctip = 'Binding species'
        self.rtbind = mQLabel('Molecule:',ctip,'Cr',12)
        self.etbind = mQLineEdit('',ctip,'l',12)
        self.rtbind.setDisabled(True)
        self.etbind.setDisabled(True)
        self.grid.addWidget(self.rtbind,8,34,1,1)
        self.grid.addWidget(self.etbind,8,35,1,2)
        # name of binding molecule from SMILES
        ctip = 'Name of binding molecule using SMILES'
        self.rtbsmi = mQLabel('Name:',ctip,'Cr',12)
        self.etbsmi = mQLineEdit('',ctip,'l',12)
        self.rtbsmi.setDisabled(True)
        self.etbsmi.setDisabled(True)
        self.grid.addWidget(self.rtbsmi,9,34,1,1)
        self.grid.addWidget(self.etbsmi,9,35,1,1)
        # separate in xyz file
        ctip = 'Separate molecules in xyz or input file with ------'
        self.chsep = mQCheckBox('separate',ctip,12)
        self.chsep.setDisabled(True)    
        self.grid.addWidget(self.chsep,9,36,1,2)
        # number of binding conformations to generate
        ctip = 'Number of different conformations to be generated'
        self.rtnbind = mQLabel('Conformations:',ctip,'Cr',12)
        self.etnbind = mQLineEdit('',ctip,'l',12)
        self.rtnbind.setDisabled(True)
        self.etnbind.setDisabled(True)
        self.grid.addWidget(self.rtnbind,10,34,1,1)
        self.grid.addWidget(self.etnbind,10,35,1,1)
        # charge of binding species
        ctip = 'Charge of binding species'
        self.rtchbind = mQLabel('Charge:',ctip,'Cr',12)
        self.etchbind = mQLineEdit('',ctip,'l',12)
        self.rtchbind.setDisabled(True)
        self.etchbind.setDisabled(True)
        self.grid.addWidget(self.rtchbind,10,36,1,1)
        self.grid.addWidget(self.etchbind,10,37,1,1)
        # min/max distance
        ctip = 'Specify placing minimum/maximum distance (in A) and/or axial/equatorial orientation'
        self.rtplace = mQLabel('Distance:',ctip,'Cr',12)
        ctip = 'Minimum distance between the two molecules. 0 corresponds the marginally non-overlapping configuration'
        self.etplacemin = mQLineEdit('',ctip,'l',12)
        ctip = 'Maximum distance between the two molecules. 0 corresponds the marginally non-overlapping configuration'
        self.etplacemax = mQLineEdit('',ctip,'l',12)
        self.rtplace.setDisabled(True)
        self.etplacemin.setDisabled(True)
        self.etplacemax.setDisabled(True)
        self.grid.addWidget(self.rtplace,11,34,1,1)
        self.grid.addWidget(self.etplacemin,11,35,1,1)
        self.grid.addWidget(self.etplacemax,11,36,1,1)
        # mask for atom/center of mass reference
        ctip = 'Reference atoms in extra molecules to be used for placement(e.g. 1,2 or 1-6 or COM or Fe) Default COM (center mass)'
        #self.rtmaskbind = mQLabel('Reference:',ctip,'r',14)
        self.etmaskbind = mQLineEdit('COM',ctip,'l',12)
        #self.rtmaskbind.setDisabled(True)
        self.etmaskbind.setDisabled(True)
        #self.grid.addWidget(self.rtmaskbind,11,37,1,1)
        self.grid.addWidget(self.etmaskbind,11,37,1,1)
        # angle/orientation
        ctip = 'Specify placement type or angle. Angle overwrites placement.'
        self.rtplacea = mQLabel('Angle:',ctip,'Cr',12)
        ctip = 'Azimouthal angle phi from 0 to 180'
        self.etplacephi = mQLineEdit('',ctip,'l',12)
        ctip = 'Polar angle theta from 0 to 360'
        self.etplacetheta = mQLineEdit('',ctip,'l',12)
        self.rtplacea.setDisabled(True)
        self.etplacephi.setDisabled(True)
        self.etplacetheta.setDisabled(True)
        self.grid.addWidget(self.rtplacea,12,34,1,1)
        self.grid.addWidget(self.etplacephi,12,35,1,1)
        self.grid.addWidget(self.etplacetheta,12,36,1,1)
        # placement of extr molecule
        ctip = 'Orientation for placing additional molecule'
        qcav = ['','axial','equatorial']
        self.dmolp = mQComboBox(qcav,ctip,12)
        self.dmolp.setDisabled(True)
        self.grid.addWidget(self.dmolp,12,37,1,1)
        # input file generation
        ctip = 'Generate input files'
        self.chkI = mQCheckBox('Input files',ctip,12)
        self.chkI.stateChanged.connect(self.enableqeinput)
        self.grid.addWidget(self.chkI,13,34,1,2)
        # jobscript generation
        ctip = 'Generate jobscripts'
        self.chkJ = mQCheckBox('Jobscripts',ctip,12)
        self.chkJ.stateChanged.connect(self.enablejinput)
        self.grid.addWidget(self.chkJ,13,36,1,2)
        # input for QC calculation
        ctip = 'Enter input for Quantum Chemistry calculations'
        self.butQc = mQPushButton('Enter QC input',ctip,12)
        self.butQc.setDisabled(True)
        self.butQc.clicked.connect(self.qcinput)
        ctip = 'Select QC code'
        qcav = ['TeraChem','GAMESS','QChem']
        self.qcode = mQComboBox(qcav,ctip,12)
        self.qcode.setDisabled(True)
        self.grid.addWidget(self.butQc,14,34,1,2)
        self.grid.addWidget(self.qcode,15,34,1,2)
        # input for jobscripts
        ctip = 'Enter input for jobscript files'
        self.butJob = mQPushButton('Enter job input',ctip,12)
        self.butJob.setDisabled(True)
        self.butJob.clicked.connect(self.jobenable)
        ctip = 'Select job scheduler'
        qcav = ['SGE','SLURM']
        self.scheduler = mQComboBox(qcav,ctip,12)
        self.scheduler.setDisabled(True)
        self.grid.addWidget(self.butJob,14,36,1,2)
        self.grid.addWidget(self.scheduler,15,36,1,2)
        # quit button
        ctip = 'Quit program'
        self.butQ = mQPushButton('Quit',ctip,14)
        self.butQ.clicked.connect(self.qexit)
        self.grid.addWidget(self.butQ,17,37,1,1)
        ################################################
        ################################################
        ################################################
        ##########################
        ### information window ###
        ##########################
        self.iWind = QWidget()
        self.iWind.setWindowTitle('Running')
        self.iWtxt = mQTextEdit('Program started..','l',14)
        self.iWtxt.setParent(self.iWind)
        self.sgrid.addWidget(self.iWind)
        #######################################
        ### create terachem-qc input window ###
        #######################################
        self.qctWindow = QWidget() # TC QC window
        self.qctgrid = QGridLayout()
        self.qctWindow.setWindowTitle('Terachem Input')
        self.sgrid.addWidget(self.qctWindow) # add to stacked grid
        self.qctWindow.setPalette(p) # set background color
        self.qctWindow.setLayout(self.qctgrid) # set layout
        c0 = mQPixmap(globs.installdir+'/icons/petachem.png')
        self.qctgrid.addWidget(c0,0,2,1,1)
        # top text
        self.txtqc = mQLabel('   TeraChem Input','','C',18)
        self.qctgrid.addWidget(self.txtqc,0,3,1,2)
        # text for specifying charge
        ctip = 'Charge of the system, default: 0'
        self.rtqctch = mQLabel('Charge:',ctip,'r',14)
        self.etqctch = mQLineEdit('0',ctip,'l',14)
        self.qctgrid.addWidget(self.rtqctch,1,1,1,1)
        self.qctgrid.addWidget(self.etqctch,1,2,1,1)
        # text for specifying spin state
        ctip = 'Spin multiplicity of the system, default: 1'
        self.rtqctspin = mQLabel('Spin:',ctip,'r',14)
        self.etqctspin = mQLineEdit('1',ctip,'l',14)
        self.qctgrid.addWidget(self.rtqctspin,1,3,1,1)
        self.qctgrid.addWidget(self.etqctspin,1,4,1,1)
        # drop menu for selecting type of calculation
        qcav = ['energy','minimize','gradient','ts']
        ctip = 'Specify calculation type, default: minimize'
        self.rtqctcalc = mQLabel('Calculation:',ctip,'r',14)
        self.qctcalc = mQComboBox(qcav,ctip,14)
        self.qctgrid.addWidget(self.rtqctcalc,2,1,1,1)
        self.qctgrid.addWidget(self.qctcalc,2,2,1,1)
        # text for specifying electronic structure method
        ctip = 'Select electronic structure method, default: ub3lyp'
        self.rtqctmethod = mQLabel('Method:',ctip,'r',14)
        self.etqctmethod = mQLineEdit('ub3lyp',ctip,'l',14)
        self.qctgrid.addWidget(self.rtqctmethod,2,3,1,1)
        self.qctgrid.addWidget(self.etqctmethod,2,4,1,1)
        # text for specifying basis set
        ctip = 'Select basis set, default: lacvp_s'
        self.rtqctbasis = mQLabel('Basis:',ctip,'r',14)
        self.etqctbasis = mQLineEdit('lacvps_ecp',ctip,'l',14)
        self.qctgrid.addWidget(self.rtqctbasis,3,1,1,1)
        self.qctgrid.addWidget(self.etqctbasis,3,2,1,1)
        # drop menu for selecting dispersion
        qcav = ['yes','no','d2','d3']
        ctip = 'Select dispersion correction'
        self.rtcdisp = mQLabel('Dispersion:',ctip,'r',14)
        self.qctsel = mQComboBox(qcav,ctip,14)
        self.qctgrid.addWidget(self.rtcdisp,3,3,1,1)
        self.qctgrid.addWidget(self.qctsel,3,4,1,1)
        self.qctsel.setCurrentIndex(1)
        # editor for additional input
        ctip='Specify additional input here'
        self.rtqctadd = mQLabel('Additional input:',ctip,'Cr',14)
        self.qceditor = mQTextEdit('','l',14)
        self.qctgrid.addWidget(self.rtqctadd,4,2,1,1)
        self.qctgrid.addWidget(self.qceditor,4,3,1,2)
        # button for addition
        ctip = 'make default'
        self.butqctlf = mQPushButton('Make default',ctip,14)
        self.butqctlf.clicked.connect(self.qctdef)
        self.qctgrid.addWidget(self.butqctlf,5,1,2,1)
        # button for addition
        ctip = 'Submit input for Quantum Chemistry'
        self.butqctSub = mQPushButton('Submit',ctip,14)
        self.butqctSub.clicked.connect(self.qretmain)
        self.qctgrid.addWidget(self.butqctSub,5,2,2,1)
        # button for return
        ctip = 'Return to main menu'
        self.butqctRet = mQPushButton('Return',ctip,14)
        self.butqctRet.clicked.connect(self.qretmain)
        self.qctgrid.addWidget(self.butqctRet,5,4,2,1)
        # load defaults if existing
        if glob.glob(globs.installdir+'/Data/.tcdefinput.inp'):
                loadfrominputtc(self,globs.installdir+'/Data/.tcdefinput.inp')
        #####################################
        ### create gamess-qc input window ###
        #######################################
        self.qcgWindow = QWidget() # TC QC window
        self.qcggrid = QGridLayout()
        self.qcgWindow.setWindowTitle('GAMESS Input')
        self.sgrid.addWidget(self.qcgWindow) # add to stacked grid
        self.qcgWindow.setPalette(p) # set background color
        self.qcgWindow.setLayout(self.qcggrid) # set layout
        # top text
        self.txtqcg = mQLabel('GAMESS Input','','c',18)
        self.qcggrid.addWidget(self.txtqcg,0,0,1,4)
        # text for specifying charge
        ctip = 'Charge of the system, default: 0'
        self.rtqcgch = mQLabel('Charge:',ctip,'r',14)
        self.etqcgch = mQLineEdit('0',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcgch,1,0,1,1)
        self.qcggrid.addWidget(self.etqcgch,1,1,1,1)
        # text for specifying spin state
        ctip = 'Spin multiplicity of the system, default: 1'
        self.rtqcgspin = mQLabel('Spin:',ctip,'r',14)
        self.etqcgspin = mQLineEdit('1',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcgspin,1,2,1,1)
        self.qcggrid.addWidget(self.etqcgspin,1,3,1,1)
        # drop menu for selecting type of calculation
        qcav = ['energy','minimize','ts']
        ctip = 'Specify calculation type, default: minimize'
        self.rtqcgcalc = mQLabel('Calculation:',ctip,'r',14)
        self.qcgcalc = mQComboBox(qcav,ctip,14)
        self.qcggrid.addWidget(self.rtqcgcalc,2,0,1,1)
        self.qcggrid.addWidget(self.qcgcalc,2,1,1,1)
        # text for specifying electronic structure method
        ctip = 'Select electronic structure method, default: ub3lyp'
        self.rtqcgmethod = mQLabel('Method:',ctip,'r',14)
        self.etqcgmethod = mQLineEdit('ub3lyp',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcgmethod,2,2,1,1)
        self.qcggrid.addWidget(self.etqcgmethod,2,3,1,1)
        # text for specifying basis set
        ctip = 'Select GBASIS input, default: 6'
        self.rtqcgbasis = mQLabel('GBASIS:',ctip,'r',14)
        self.etqcgbasis = mQLineEdit('6',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcgbasis,3,0,1,1)
        self.qcggrid.addWidget(self.etqcgbasis,3,1,1,1)
        # text for specifying basis set
        ctip = 'Select NGAUSS input, default: N31'
        self.rtqcngauss = mQLabel('NGAUSS:',ctip,'r',14)
        self.etqcngauss = mQLineEdit('N31',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcngauss,3,2,1,1)
        self.qcggrid.addWidget(self.etqcngauss,3,3,1,1)
        # text for specifying polarization functions
        ctip = 'Select NPFUNC input, default: 0'
        self.rtqcnpfunc = mQLabel('NPFUNC:',ctip,'r',14)
        self.etqcnpfunc = mQLineEdit('',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcnpfunc,4,0,1,1)
        self.qcggrid.addWidget(self.etqcnpfunc,4,1,1,1)
        # text for specifying polarization functions
        ctip = 'Select NDFUNC input, default: 0'
        self.rtqcndfunc = mQLabel('NDFUNC:',ctip,'r',14)
        self.etqcndfunc = mQLineEdit('',ctip,'l',14)
        self.qcggrid.addWidget(self.rtqcndfunc,4,2,1,1)
        self.qcggrid.addWidget(self.etqcndfunc,4,3,1,1)
        # editor for additional input
        ctip='Specify additional input for block SYS here'
        self.rtqcgadd1 = mQLabel('SYS input:',ctip,'Cr',14)
        self.qcgedsys = mQTextEdit('','l',12)
        self.qcggrid.addWidget(self.rtqcgadd1,5,0,1,1)
        self.qcggrid.addWidget(self.qcgedsys,5,1,1,1)
        # editor for additional input
        ctip='Specify additional input for block CTRL here'
        self.rtqcgadd1 = mQLabel('CTRL input:',ctip,'Cr',14)
        self.qcgedctrl = mQTextEdit('','l',12)
        self.qcggrid.addWidget(self.rtqcgadd1,5,2,1,1)
        self.qcggrid.addWidget(self.qcgedctrl,5,3,1,1)
        # editor for additional input
        ctip='Specify additional input for block SCF here'
        self.rtqcgadd3 = mQLabel('SCF input:',ctip,'Cr',14)
        self.qcgedscf = mQTextEdit('','l',12)
        self.qcggrid.addWidget(self.rtqcgadd3,6,0,1,1)
        self.qcggrid.addWidget(self.qcgedscf,6,1,1,1)
        # editor for additional input
        ctip='Specify additional input for block STAT here'
        self.rtqcgadd4 = mQLabel('STAT input:',ctip,'Cr',14)
        self.qcgedstat = mQTextEdit('','l',12)
        self.qcggrid.addWidget(self.rtqcgadd4,6,2,1,1)
        self.qcggrid.addWidget(self.qcgedstat,6,3,1,1)
        # button for addition
        ctip = 'make default'
        self.butqcglf = mQPushButton('Make default',ctip,14)
        self.butqcglf.clicked.connect(self.qcgdef)
        self.qcggrid.addWidget(self.butqcglf,7,0,1,1)
        # button for addition
        ctip = 'Submit input for Quantum Chemistry'
        self.butqcgSub = mQPushButton('Submit',ctip,14)
        self.butqcgSub.clicked.connect(self.qretmain)
        self.qcggrid.addWidget(self.butqcgSub,7,1,1,2)
        # button for return
        ctip = 'Return to main menu'
        self.butqcgRet = mQPushButton('Return',ctip,14)
        self.butqcgRet.clicked.connect(self.qretmain)
        self.qcggrid.addWidget(self.butqcgRet,7,3,1,1)
        self.qcggrid.setRowMinimumHeight(7,30)
        # load defaults if existing
        if glob.glob(globs.installdir+'/Data/.gamdefinput.inp'):
                loadfrominputgam(self,globs.installdir+'/Data/.gamdefinput.inp')
        #######################################
        #### create Qchem-qc input window #####
        #######################################
        self.qcQWindow = QWidget() # TC QC window
        self.qcQgrid = QGridLayout()
        self.qctWindow.setWindowTitle('QChem Input')
        self.sgrid.addWidget(self.qcQWindow) # add to stacked grid
        self.qcQWindow.setPalette(p) # set background color
        self.qcQWindow.setLayout(self.qcQgrid) # set layout
        # top text
        self.txtQqc = mQLabel('QChem Input','','c',18)
        self.qcQgrid.addWidget(self.txtQqc,0,0,1,5)
        # text for specifying charge
        ctip = 'Charge of the system, default: 0'
        self.rtqcQch = mQLabel('Charge:',ctip,'r',14)
        self.etqcQch = mQLineEdit('0',ctip,'l',14)
        self.qcQgrid.addWidget(self.rtqcQch,1,0,1,1)
        self.qcQgrid.addWidget(self.etqcQch,1,1,1,1)
        # text for specifying spin state
        ctip = 'Spin multiplicity of the system, default: 1'
        self.rtqcQspin = mQLabel('Spin:',ctip,'r',14)
        self.etqcQspin = mQLineEdit('1',ctip,'l',14)
        self.qcQgrid.addWidget(self.rtqcQspin,1,2,1,1)
        self.qcQgrid.addWidget(self.etqcQspin,1,3,1,1)
        # additional molecule
        ctip = 'Unrestricted calculation?'
        self.chQun = mQCheckBox('Unrestricted',ctip,14)
        self.chQun.setChecked(True)
        self.qcQgrid.addWidget(self.chQun,1,4,1,1)
        # drop menu for selecting type of calculation
        qcav = ['energy','minimize','ts']
        ctip = 'Specify calculation type, default: minimize'
        self.rtqcQcalc = mQLabel('Calculation:',ctip,'r',14)
        self.qcQcalc = mQComboBox(qcav,ctip,14)
        self.qcQgrid.addWidget(self.rtqcQcalc,2,0,1,1)
        self.qcQgrid.addWidget(self.qcQcalc,2,1,1,1)
        # text for specifying basis set
        ctip = 'Select basis set, default: lanl2dz'
        self.rtqcQbasis = mQLabel('Basis:',ctip,'r',14)
        self.etqcQbasis = mQLineEdit('lanl2dz',ctip,'l',14)
        self.qcQgrid.addWidget(self.rtqcQbasis,2,2,1,1)
        self.qcQgrid.addWidget(self.etqcQbasis,2,3,1,1)
        # text for specifying exchange
        ctip = 'Select exchange, default: b3lyp'
        self.rtqcQex = mQLabel('Exchange:',ctip,'r',14)
        self.etqcQex = mQLineEdit('b3lyp',ctip,'l',14)
        self.qcQgrid.addWidget(self.rtqcQex,3,0,1,1)
        self.qcQgrid.addWidget(self.etqcQex,3,1,1,1)
        # text for specifying electronic structure method
        ctip = 'Select correlation, default: none'
        self.rtqcQcor = mQLabel('Correlation:',ctip,'r',14)
        self.etqcQcor = mQLineEdit('none',ctip,'l',14)
        self.qcQgrid.addWidget(self.rtqcQcor,3,2,1,1)
        self.qcQgrid.addWidget(self.etqcQcor,3,3,1,1)
        # editor for additional input
        ctip='Specify additional input here'
        self.rtqcQadd = mQLabel('Additional input:',ctip,'Cr',14)
        self.qcQeditor = mQTextEdit('','l',12)
        self.qcQgrid.addWidget(self.rtqcQadd,4,1,1,1)
        self.qcQgrid.addWidget(self.qcQeditor,4,2,1,2)
        # button for addition
        ctip = 'make default'
        self.butqcQlf = mQPushButton('Make default',ctip,14)
        self.butqcQlf.clicked.connect(self.qcqdef)
        self.qcQgrid.addWidget(self.butqcQlf,5,0,1,1)
        # button for addition
        ctip = 'Submit input for Quantum Chemistry'
        self.butqcQSub = mQPushButton('Submit',ctip,14)
        self.butqcQSub.clicked.connect(self.qretmain)
        self.qcQgrid.addWidget(self.butqcQSub,5,2,1,2)
        # button for return
        ctip = 'Return to main menu'
        self.butqcQRet = mQPushButton('Return',ctip,14)
        self.butqcQRet.clicked.connect(self.qretmain)
        self.qcQgrid.addWidget(self.butqcQRet,5,4,1,1)
        # load defaults if existing
        if glob.glob(globs.installdir+'/Data/.qchdefinput.inp'):
            loadfrominputqch(self,globs.installdir+'/Data/.qchdefinput.inp')
        #####################################
        ### create jobscript input window ###
        #####################################
        self.jWindow = QWidget() # jobscript window
        self.jgrid = QGridLayout()
        self.jWindow.setWindowTitle('Jobscript parameters')
        self.sgrid.addWidget(self.jWindow) # add to stacked grid
        self.jWindow.setPalette(p) # set background color
        self.jWindow.setLayout(self.jgrid) # set layout
        c1 = mQPixmap(globs.installdir+'/icons/sge.png')
        c2 = mQPixmap(globs.installdir+'/icons/slurm.png')
        self.jgrid.addWidget(c1,1,2,1,1)
        self.jgrid.addWidget(c2,1,4,1,1)
        # top text
        self.txtj = mQLabel('Jobscript Input','','r',18)
        self.jgrid.addWidget(self.txtj,0,2,1,2)
        # text for main job identifier
        ctip = 'Job main identifier, e.g. feII'
        self.rtjname = mQLabel('Job name:',ctip,'r',14)
        self.etjname = mQLineEdit('myjob',ctip,'l',14)
        self.jgrid.addWidget(self.rtjname,2,1,1,1)
        self.jgrid.addWidget(self.etjname,2,2,1,1)
        # text for specifying spin state
        ctip = 'Queue to use, e.g. gpus'
        self.rtjqueue = mQLabel('Queue:',ctip,'r',14)
        self.etjqueue = mQLineEdit('',ctip,'l',14)
        self.jgrid.addWidget(self.rtjqueue,2,3,1,1)
        self.jgrid.addWidget(self.etjqueue,2,4,1,1)
        # text for specifying wall time
        ctip = 'Wall time request, e.g. 48'
        self.rtjwallt = mQLabel('Wall time:',ctip,'r',14)
        self.etjwallt = mQLineEdit('48h',ctip,'l',14)
        self.jgrid.addWidget(self.rtjwallt,3,1,1,1)
        self.jgrid.addWidget(self.etjwallt,3,2,1,1)
        # text for specifying memory
        ctip = 'Memory request, e.g. 10G'
        self.rtjmem = mQLabel('Memory:',ctip,'r',14)
        self.etjmem = mQLineEdit('10G',ctip,'l',14)
        self.jgrid.addWidget(self.rtjmem,3,3,1,1)
        self.jgrid.addWidget(self.etjmem,3,4,1,1)
        # text for specifying charge
        ctip = 'Number of CPUs requested, default: 0'
        self.rtjcpus = mQLabel('CPUs:',ctip,'r',14)
        self.etjcpus = mQLineEdit('',ctip,'l',14)
        self.jgrid.addWidget(self.rtjcpus,4,1,1,1)
        self.jgrid.addWidget(self.etjcpus,4,2,1,1)
        # text for specifying spin state
        ctip = 'Number of GPUs requested, default: 0'
        self.rtjgpus = mQLabel('GPUs:',ctip,'r',14)
        self.etjgpus = mQLineEdit('',ctip,'l',14)
        self.jgrid.addWidget(self.rtjgpus,4,3,1,1)
        self.jgrid.addWidget(self.etjgpus,4,4,1,1)
        # text for modules to be loaded
        ctip = 'Modules to be loaded, e.g. terachem, openmpi'
        self.rtjmod = mQLabel('Modules:',ctip,'r',14)
        self.etjmod = mQLineEdit('',ctip,'l',14)
        self.jgrid.addWidget(self.rtjmod,5,1,1,1)
        self.jgrid.addWidget(self.etjmod,5,2,1,1)
        # editor for additional input
        ctip='Specify additional input for initial options, e.g. -j y'
        self.rtjadd1 = mQLabel('Options:',ctip,'r',14)
        self.etjopt = mQTextEdit('','l',12)
        self.jgrid.addWidget(self.rtjadd1,6,1,1,1)
        self.jgrid.addWidget(self.etjopt,6,2,1,1)
        # editor for additional input
        ctip='Specify additional commands, e.g. mkdir $WORKDIR/test'
        self.rtjadd2 = mQLabel('Commands:',ctip,'Cr',14)
        self.jcomm = mQTextEdit('','l',12)
        self.jgrid.addWidget(self.rtjadd2,6,3,1,1)
        self.jgrid.addWidget(self.jcomm,6,4,1,1)
        # button for addition
        ctip = 'make default'
        self.butqcJlf = mQPushButton('Make default',ctip,14)
        self.butqcJlf.clicked.connect(self.jobdef)
        self.jgrid.addWidget(self.butqcJlf,7,1,1,1)
        # button for addition
        ctip = 'Submit input for Quantum Chemistry'
        self.butqcgSub = mQPushButton('Submit',ctip,14)
        self.butqcgSub.clicked.connect(self.qretmain)
        self.jgrid.addWidget(self.butqcgSub,7,2,1,2)
        # button for return
        ctip = 'Return to main menu'
        self.butqcgRet = mQPushButton('Return',ctip,14)
        self.butqcgRet.clicked.connect(self.qretmain)
        self.jgrid.addWidget(self.butqcgRet,7,4,1,1)
        # load defaults if existing
        if glob.glob(globs.installdir+'/Data/.jobdefinput.inp'):
            loadfrominputjob(self,globs.installdir+'/Data/.jobdefinput.inp')
        ##########################################
        ### create local DB interaction window ###
        ##########################################
        self.DBWindow = QWidget() # DB Window
        self.DBlgrid = QGridLayout()
        self.sgrid.addWidget(self.DBWindow)
        self.DBWindow.setPalette(p)
        self.DBWindow.setLayout(self.DBlgrid)
        self.DBWindow.setWindowTitle('Insert/remove to/from Database')
        # top text
        self.txtdb = mQLabel('Database Update','','c',18)
        self.DBlgrid.addWidget(self.txtdb,0,0,1,4)
        # drop menu for selecting type
        ctip = 'Select what type of molecule you want to add/remove to/from the database'
        self.rtDBsel = mQLabel('Select type:','','r',14)
        qcav = ['core','ligand','binding']
        self.DBsel = mQComboBox(qcav,ctip,14)
        self.DBsel.setCurrentIndex(1)
        self.DBlgrid.addWidget(self.rtDBsel,1,1,1,1)
        self.DBlgrid.addWidget(self.DBsel,1,2,1,1)
        # text for selecting type
        ctip = 'Type SMILES string for molecule'
        self.rtDBsmi = mQLabel('SMILES or file:','','r',14)
        self.etDBsmi = mQLineEdit('',ctip,'l',14)
        self.DBlgrid.addWidget(self.rtDBsmi,2,1,1,1)
        self.DBlgrid.addWidget(self.etDBsmi,2,2,1,1)
        # text for specifying name
        ctip = 'Type name for molecule'
        self.rtDBname = mQLabel('Name:','','r',14)
        self.etDBname = mQLineEdit('',ctip,'l',14)
        self.DBlgrid.addWidget(self.rtDBname,3,1,1,1)
        self.DBlgrid.addWidget(self.etDBname,3,2,1,1)
        # text for specifying category
        ctip = 'Type groups for ligand'
        self.rtDBgrps = mQLabel('Groups:','','r',14)
        self.etDBgrps = mQLineEdit('',ctip,'l',14)
        self.DBlgrid.addWidget(self.rtDBgrps,4,1,1,1)
        self.DBlgrid.addWidget(self.etDBgrps,4,2,1,1)
        # text for specifying groups
        ctip = 'Type category for ligand. It can be used for building complexes("build")\n'
        ctip += 'functionalizing existing cores("functionalize") or both.'
        self.rtDBctg = mQLabel('Category:','','r',14)
        self.etDBctg = mQComboBox(['all','build','functionalize'],ctip,14)
        self.DBlgrid.addWidget(self.rtDBctg,5,1,1,1)
        self.DBlgrid.addWidget(self.etDBctg,5,2,1,1)
        # checkboxes for FF optimization
        ctip = 'Force-field optimize in isolation'
        self.lFFb = mQCheckBox('FF before',ctip,14)
        self.DBlgrid.addWidget(self.lFFb,6,1,1,1)
        self.lFFb.setChecked(True)
        # checkboxes for FF optimization
        ctip = 'Force-field optimize in molecule'
        self.lFFa = mQCheckBox('FF after',ctip,14)
        self.DBlgrid.addWidget(self.lFFa,6,2,1,1)
        self.lFFa.setChecked(True)
        # drop menu for denticity
        ctip = 'Type denticity for SMILES molecule'
        self.rtDBsmident = mQLabel('Denticity:','','r',14)
        qcav = ['1','2','3','4','5','6','7','8','9','10']
        self.DBdent = mQComboBox(qcav,ctip,14)
        self.DBdent.setCurrentIndex(0)
        self.DBlgrid.addWidget(self.rtDBsmident,7,1,1,1)
        self.DBlgrid.addWidget(self.DBdent,7,2,1,1)
        self.DBsel.currentIndexChanged.connect(self.dbchange)  ### error if not here
        # text for typing input connection atoms 
        ctip = 'Type indices for connection atoms, default: 1'
        self.rtDBsmicat = mQLabel('Catoms:',ctip,'r',14)
        self.etDBsmicat = mQLineEdit('',ctip,'l',14)
        self.DBlgrid.addWidget(self.rtDBsmicat,8,1,1,1)
        self.DBlgrid.addWidget(self.etDBsmicat,8,2,1,1)
        # button for addition
        ctip = 'Load molecule from file'
        self.butDBAlf = mQPushButton('Load file..',ctip,14)
        self.butDBAlf.clicked.connect(self.qDBload)
        self.DBlgrid.addWidget(self.butDBAlf,9,0,1,1)
        # button for addition
        ctip = 'Add new molecule to database'
        self.butDBAub = mQPushButton('Add',ctip,14)
        self.butDBAub.clicked.connect(self.qaddDB)
        self.DBlgrid.addWidget(self.butDBAub,9,1,1,1)
        # button for removal
        ctip = 'Remove molecule from database'
        self.butDBDub = mQPushButton('Remove',ctip,14)
        self.butDBDub.clicked.connect(self.qdelDB)
        self.DBlgrid.addWidget(self.butDBDub,9,2,1,1)
        # button for return
        ctip = 'Return to main menu'
        self.butDBRet = mQPushButton('Return',ctip,14)
        self.butDBRet.clicked.connect(self.qretmain)
        self.DBlgrid.addWidget(self.butDBRet,9,3,1,1)
        #############################################
        ### create chemical DB interaction window ###
        #############################################
        self.cDBWindow = QWidget() # DB Window
        self.cDBgrid = QGridLayout()
        self.sgrid.addWidget(self.cDBWindow)
        self.cDBWindow.setPalette(p)
        self.cDBWindow.setLayout(self.cDBgrid)
        self.cDBWindow.setWindowTitle('Chemical Database Search')
        ctip = 'Specify screening options or similarity search.'
        # top text
        self.txtcdb = mQLabel('Database Search',ctip,'c',18)
        self.cDBgrid.addWidget(self.txtcdb,0,0,1,9)
        # text for reference
        self.rtcDBref = mQLabel('Similarity search','','c',18)
        self.cDBgrid.addWidget(self.rtcDBref,1,0,1,3)
        # text for selecting type
        ctip = 'Type SMILES string or molecule name for reference molecule in similarity search'
        self.rtcDBsmi = mQLabel('SMILES:',ctip,'r',14)
        self.etcDBsmi = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.rtcDBsmi,2,0,1,1)
        self.cDBgrid.addWidget(self.etcDBsmi,2,1,1,2)
        # button for loading from file
        ctip = 'Load reference molecule from file'
        self.rtcDBAlf = mQLabel('From file:',ctip,'r',14)
        self.butcDBAlf = mQPushButton('Load..',ctip,14)
        self.butcDBAlf.clicked.connect(self.qcDBload)
        self.cDBgrid.addWidget(self.rtcDBAlf,3,0,1,1)
        self.cDBgrid.addWidget(self.butcDBAlf,3,1,1,2)
        # connection atoms for smarts/smiles
        ctip = 'Specify the connection atoms in SMARTS/SMILES. Default: 1'
        self.rtcDBcatoms = mQLabel('Conn atoms:',ctip,'r',14)
        self.etcDBcatoms = mQLineEdit('1',ctip,'l',14)
        self.cDBgrid.addWidget(self.rtcDBcatoms,4,0,1,1)
        self.cDBgrid.addWidget(self.etcDBcatoms,4,1,1,1)
        # how many molecules to return
        ctip = 'Specify the number of results you want.'
        self.rtcDBres = mQLabel('Results:',ctip,'r',14)
        self.etcDBnres = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.rtcDBres,5,0,1,1)
        self.cDBgrid.addWidget(self.etcDBnres,5,1,1,1)
        # text for output file
        ctip = 'Please type in output file'
        self.rtcDBsoutf = mQLabel('Output file:',ctip,'r',14)
        self.etcDBoutf = mQLineEdit('simres',ctip,'l',14)
        self.cDBgrid.addWidget(self.rtcDBsoutf,6,0,1,1)
        self.cDBgrid.addWidget(self.etcDBoutf,6,1,1,1)
        # drop menu for output file
        qcav = ['.smi']#,'.mol','.sdf']
        self.cDBdent = mQComboBox(qcav,ctip,14)
        self.cDBgrid.addWidget(self.cDBdent,6,2,1,1)
        # text for screening options
        self.rtcDBsc = mQLabel('Screening options',ctip,'c',18)
        self.cDBgrid.addWidget(self.rtcDBsc,1,6,1,4)
        # text for selecting database
        self.rtcDBsel = mQLabel('Select database:',ctip,'c',14)
        self.cDBgrid.addWidget(self.rtcDBsel,2,5,1,2)
        # get existing databases
        ctip = 'Select the database you want to use. Please note that fastsearch indexes DBs allow for much faster screening.'
        dbdir = globs.chemdbdir
        dbs0 = glob.glob(dbdir+"/*.sdf")
        dbs1 = [d.rsplit('/',1)[-1] for d in dbs0]
        dbs = [d.split('.',1)[0] for d in dbs1]
        self.cDBsel = mQComboBox(dbs,ctip,14)
        self.cDBgrid.addWidget(self.cDBsel,2,7,1,2)
        # text for SMARTS pattern
        ctip = 'Type SMARTS pattern for matching'
        self.rtcDBsmarts = mQLabel('SMARTS:',ctip,'r',14)
        self.etcDBsmarts = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.rtcDBsmarts,3,5,1,2)
        self.cDBgrid.addWidget(self.etcDBsmarts,3,7,1,2)
        ctip = 'Select molecular fingerprint'
        # text for selecting fingerprint
        self.rtcDBsf = mQLabel('Select Fingerprint:',ctip,'r',14)
        self.cDBgrid.addWidget(self.rtcDBsf,4,5,1,2)
        # select fingerprint
        opts = ['FP2','FP3','FP4','MACCS']
        self.cDBsf = mQComboBox(opts,ctip,14)
        self.cDBgrid.addWidget(self.cDBsf,4,7,1,2)
        # get options for screening
        ctip = 'Specify minimum and maximum values for filters.'
        self.rtcDBmin = mQLabel('min',ctip,'c',14)
        self.rtcDBmax = mQLabel('max',ctip,'c',14)
        self.cDBgrid.addWidget(self.rtcDBmin,6,7,1,1)
        self.cDBgrid.addWidget(self.rtcDBmax,6,8,1,1)
        ctip = 'Total number of atoms.'
        self.rtcDBat0 = mQLabel('atoms:',ctip,'r',14)
        self.cDBgrid.addWidget(self.rtcDBat0,7,5,1,2)
        ctip = 'Minimum number of atoms.'
        self.etcDBsatoms0 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsatoms0,7,7,1,1)
        ctip = 'Maximum number of atoms.'
        self.etcDBsatoms1 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsatoms1,7,8,1,1)
        ctip = 'Total number of bonds.'
        self.rtcDBb0 = mQLabel('bonds:',ctip,'r',14)
        self.cDBgrid.addWidget(self.rtcDBb0,8,5,1,2)
        ctip = 'Minimum number of total bonds.'
        self.etcDBsbonds0 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsbonds0,8,7,1,1)
        ctip = 'Maximum number of total bonds.'
        self.etcDBsbonds1 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsbonds1,8,8,1,1)
        ctip = 'Total number of aromatic.'
        self.rtcDBba0 = mQLabel('aromatic bonds:',ctip,'r',14)
        self.cDBgrid.addWidget(self.rtcDBba0,9,5,1,2)
        ctip = 'Minimum number of aromatic bonds.'
        self.etcDBsabonds0 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsabonds0,9,7,1,1)
        ctip = 'Maximum number of aromatic bonds.'
        self.etcDBsabonds1 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsabonds1,9,8,1,1)
        ctip = 'Total number of single bonds.'
        self.rtcDBs0 = mQLabel('single bonds:',ctip,'r',14)
        self.cDBgrid.addWidget(self.rtcDBs0,10,5,1,2)
        ctip = 'Minimum number of single bonds.'
        self.etcDBsbondss0 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsbondss0,10,7,1,1)
        ctip = 'Maximum number of single bonds.'
        self.etcDBsbondss1 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBsbondss1,10,8,1,1)
        ctip = 'Molecular weight.'
        self.rtcDBbtmw0 = mQLabel('MW:',ctip,'',14)
        self.cDBgrid.addWidget(self.rtcDBbtmw0,11,5,1,2)
        ctip = 'Minimum molecular weight.'
        self.etcDBmw0 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBmw0,11,7,1,1)
        ctip = 'Maximum molecular weight.'
        self.etcDBmw1 = mQLineEdit('',ctip,'l',14)
        self.cDBgrid.addWidget(self.etcDBmw1,11,8,1,1)
        # aspirin icon
        c = mQPixmap(globs.installdir+'/icons/chemdb.png')
        relresize(c,c,0.4)
        self.cDBgrid.addWidget(c,7,0,4,2)
        # button for addition
        ctip = 'Search database:'
        self.butcDBAub = mQPushButton('Search',ctip,14)
        self.butcDBAub.clicked.connect(self.qaddcDB)
        self.cDBgrid.addWidget(self.butcDBAub,11,0,1,2)
        # button for addition
        ctip = 'Draw results'
        self.butcDBd0 = mQPushButton('Draw',ctip,14)
        self.butcDBd0.clicked.connect(self.drawres)
        self.cDBgrid.addWidget(self.butcDBd0,9,2,1,1)
        # button for return
        ctip = 'Return to main menu'
        self.butcDBRet = mQPushButton('Return',ctip,14)
        self.butcDBRet.clicked.connect(self.qretmain)
        self.cDBgrid.addWidget(self.butcDBRet,11,2,1,1)
        ###########################
        ### post-process window ###
        ###########################
        self.pWindow = QWidget() # Post Window
        self.pgrid = QGridLayout()
        self.sgrid.addWidget(self.pWindow)
        self.pWindow.setPalette(p)
        self.pWindow.setLayout(self.pgrid)
        self.pWindow.setWindowTitle('Post-processing')
        self.pgrid.setRowMinimumHeight(0,50)
        self.pgrid.setRowMinimumHeight(2,30)
        # top text
        self.txtp = mQLabel('Post-processing setup','','c',18)
        self.pgrid.addWidget(self.txtp,0,0,1,4)
        # input for jobscripts
        ctip = 'Select top jobs directory'
        self.butpd = mQPushButton('Select directory',ctip,16)
        self.butpd.clicked.connect(self.pdload)
        self.pgrid.addWidget(self.butpd,1,0,1,1)
        self.etpdir = mQLineEdit(globs.rundir,ctip,'l',16)
        self.pgrid.addWidget(self.etpdir,1,1,1,1)
        # select qc code
        ctip = 'Select QC code'
        qcav = ['TeraChem','GAMESS']
        self.rpcode = mQLabel('QC code:','','c',16)
        self.pgrid.addWidget(self.rpcode,1,2,1,1)
        self.pqcode = mQComboBox(qcav,ctip,14)
        self.pgrid.addWidget(self.pqcode,1,3,1,1)
        # general summary
        ctip = 'Generate results summary'
        self.psum = mQCheckBox('Summary',ctip,16)
        self.pgrid.addWidget(self.psum,3,1,1,1)
        # metal charges
        ctip = 'Calculate metal charge'
        self.pch = mQCheckBox('Charge',ctip,16)
        self.pgrid.addWidget(self.pch,3,2,1,1)
        # wavefunction properties
        ctip = 'Generate average properties of the wavefunction'
        self.pwfnav = mQCheckBox('Wavefunction',ctip,16)
        self.pgrid.addWidget(self.pwfnav,4,1,1,1)
        # cube files
        ctip = 'Generate cubefiles'
        self.pcub = mQCheckBox('Cubes',ctip,16)
        self.pgrid.addWidget(self.pcub,4,2,1,1)
        # molecular orbital information
        ctip = 'Molecular Orbital information'
        self.porbs = mQCheckBox('MO',ctip,16)
        self.pgrid.addWidget(self.porbs,5,1,1,1)
        # NBO analysis
        ctip = 'NBO analysis'
        self.pnbo = mQCheckBox('NBO',ctip,16)
        self.pgrid.addWidget(self.pnbo,5,2,1,1)
        # d-orbital information
        #ctip = 'd-orbital information'
        #self.pdorbs = mCheck(self.pWindow,0.55,0.525,0.2,0.1,'d-orbitals',ctip,16)
        # delocalization indices
        ctip = 'Localization and delocalization indices'
        self.pdeloc = mQCheckBox('Deloc',ctip,16)
        self.pgrid.addWidget(self.pdeloc,6,1,1,1)
        # button for addition
        ctip = 'Run post-processing'
        self.butpR = mQPushButton('Run',ctip,16)
        self.pgrid.addWidget(self.butpR,7,1,1,2)
        self.butpR.clicked.connect(self.postprocGUI)
        # button for return
        ctip = 'Return to main menu'
        self.butpret = mQPushButton('Return',ctip,16)
        self.butpret.clicked.connect(self.qretmain)
        self.pgrid.addWidget(self.butpret,7,3,1,2)
        #c1p = mPic(self.pWindow,globs.installdir+'/icons/wft1.png',0.04,0.7,0.2)
        c3p = mQPixmap(globs.installdir+'/icons/wft3.png')
        self.pgrid.addWidget(c3p,3,0,4,1)
        #c2p = mPic(self.pWindow,globs.installdir+'/icons/wft2.png',0.04,0.035,0.2)
        ##################################
        ### create add geometry window ###
        ##################################
        self.geWindow = QWidget() # Geometry window
        self.gegrid = QGridLayout()
        self.sgrid.addWidget(self.geWindow)
        self.geWindow.setPalette(p)
        self.geWindow.setLayout(self.gegrid)
        self.geWindow.setWindowTitle('Add new geometry')
        # top text
        self.txtg = mQLabel('New geometry','','c',18)
        self.gegrid.addWidget(self.txtg,0,0,1,4)
        # text for specifying name
        ctip = 'Type name of geometry'
        self.rtgname = mQLabel('Name:','','r',14)
        self.etgname = mQLineEdit('',ctip,'l',14)
        self.gegrid.addWidget(self.rtgname,1,0,1,1)
        self.gegrid.addWidget(self.etgname,1,1,1,1)
        # text for specifying short name
        ctip = 'Type short identifier of geometry (2-4 letters)'
        self.rtgshort = mQLabel('Short identifier:','','r',14)
        self.etgshort = mQLineEdit('',ctip,'l',14)
        self.gegrid.addWidget(self.rtgshort,2,0,1,1)
        self.gegrid.addWidget(self.etgshort,2,1,1,1)
        # xyz file
        ctip = 'Load xyz file with geometry.'
        self.rtgf = mQLabel('Filename:','','r',14)
        self.etgf = mQLineEdit('',ctip,'l',14)
        self.gegrid.addWidget(self.rtgf,3,0,1,1)
        self.gegrid.addWidget(self.etgf,3,1,1,1)
        # load file
        self.butge = mQPushButton('Load file..',ctip,14)
        self.butge.clicked.connect(self.qgeomload)
        self.gegrid.addWidget(self.butge,4,0,1,1)
        # button for addition
        ctip = 'Add new geometry to database'
        self.butgAub = mQPushButton('Add',ctip,14)
        self.butgAub.clicked.connect(self.qaddg)
        self.gegrid.addWidget(self.butgAub,4,1,1,1)
        # button for removal
        ctip = 'Remove geometry from database'
        self.butgDub = mQPushButton('Remove',ctip,14)
        self.butgDub.clicked.connect(self.qdelg)
        self.gegrid.addWidget(self.butgDub,5,0,1,1)
        # button for return
        ctip = 'Return to main menu'
        self.butgRet = mQPushButton('Return',ctip,14)
        self.butgRet.clicked.connect(self.qretmain)
        self.gegrid.addWidget(self.butgRet,5,1,1,1)
        ################################
        ### create ligands 2D window ###
        ################################
        self.lwindow = QWidget()
        self.lwindow.setPalette(p)
        self.lgrid = QGridLayout()
        self.lwindow.setLayout(self.lgrid)
        self.sgrid.addWidget(self.lwindow)
        self.lwindow.setWindowTitle('Ligands 2D')
        self.c1p = QWidget()
        ####################
        ####################
        ####################
        ### Run main GUI ###
        ####################
        ####################
        ####################
        self.wmain.setLayout(self.grid)
        self.wmwindow.setCentralWidget(self.wmain)
        self.sgrid.setCurrentWidget(self.wmain)
        self.wmain.showMaximized()
        self.wmwindow.showMaximized()
        # resize other windows
        relresize(self.iWind,self.wmwindow,0.7)
        relresize(self.iWtxt,self.iWind,1.0)
        app.processEvents()
        app.exec_()
    '''
    #############################
    ### Callbacks for buttons ###
    #############################
    '''
    ###################################
    ###### Add new ligands input ######
    ###################################
    def addlig0(self):
        txt0 = self.lig0.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig0add.setDisabled(True)
            self.lig0add.hide()
            self.lig1add.setDisabled(False)
            self.lig1add.show()
            self.lig1.setDisabled(False)
            self.lig1occ.setDisabled(False)
            self.lig1conn.setDisabled(False)
            self.lig1H.setDisabled(False)
            self.lig1ML.setDisabled(False)
            self.lig1an.setDisabled(False)
            self.lig1nam.setDisabled(False)
    def addlig1(self):
        txt0 = self.lig1.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig1add.setDisabled(True)
            self.lig1add.hide()
            self.lig2add.setDisabled(False)
            self.lig2add.show()
            self.lig2.setDisabled(False)
            self.lig2occ.setDisabled(False)
            self.lig2conn.setDisabled(False)
            self.lig2H.setDisabled(False)
            self.lig2ML.setDisabled(False)
            self.lig2an.setDisabled(False)
            self.lig2nam.setDisabled(False)
    def addlig2(self):
        txt0 = self.lig2.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig2add.setDisabled(True)
            self.lig2add.hide()
            self.lig3add.setDisabled(False)
            self.lig3add.show()
            self.lig3.setDisabled(False)
            self.lig3occ.setDisabled(False)
            self.lig3conn.setDisabled(False)
            self.lig3H.setDisabled(False)
            self.lig3ML.setDisabled(False)
            self.lig3an.setDisabled(False)
            self.lig3nam.setDisabled(False)
    def addlig3(self):
        txt0 = self.lig3.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig3add.setDisabled(True)
            self.lig3add.hide()
            self.lig4add.setDisabled(False)
            self.lig4add.show()
            self.lig4.setDisabled(False)
            self.lig4occ.setDisabled(False)
            self.lig4conn.setDisabled(False)
            self.lig4H.setDisabled(False)
            self.lig4ML.setDisabled(False)
            self.lig4an.setDisabled(False)
            self.lig4nam.setDisabled(False)
    def addlig4(self):
        txt0 = self.lig4.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig4add.setDisabled(True)
            self.lig4add.hide()
            self.lig5add.setDisabled(False)
            self.lig5add.show()
            self.lig5.setDisabled(False)
            self.lig5occ.setDisabled(False)
            self.lig5conn.setDisabled(False)
            self.lig5H.setDisabled(False)
            self.lig5ML.setDisabled(False)
            self.lig5an.setDisabled(False)
            self.lig5nam.setDisabled(False)
    def addlig5(self):
        txt0 = self.lig0.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig5add.setDisabled(True)
            self.lig5add.hide()
            self.lig6add.setDisabled(False)
            self.lig6add.show()
            self.lig6.setDisabled(False)
            self.lig6occ.setDisabled(False)
            self.lig6conn.setDisabled(False)
            self.lig6H.setDisabled(False)
            self.lig6ML.setDisabled(False)
            self.lig6an.setDisabled(False)
            self.lig6nam.setDisabled(False)
    def addlig6(self):
        txt0 = self.lig0.text().replace(' ','')
        if len(txt0)==0:
            mQDialogWarn('No ligand specified','Please specify a ligand before adding another one.')
        else:
            self.lig6add.setDisabled(True)
            self.lig6add.hide()
            self.lig7.setDisabled(False)
            self.lig7occ.setDisabled(False)
            self.lig7conn.setDisabled(False)
            self.lig7H.setDisabled(False)
            self.lig7ML.setDisabled(False)
            self.lig7an.setDisabled(False)
            self.lig7nam.setDisabled(False)
    ##############################
    ### Local Database window ####
    ##############################
    ### load molecule from file
    def qDBload(self):
        name = QFileDialog.getOpenFileName(self.DBWindow,'Open File','.',"Molecule files *.xyz, *.mol (*.xyz *.mol)")
        if name[0] != '':
            self.etDBsmi.setText(os.path.relpath(name[0]))
    ### enable add to database interface
    def enableDB(self):
        self.DBWindow.setWindowModality(2)
        self.DBWindow.show()
        self.sgrid.setCurrentWidget(self.DBWindow)
    ### callback for addition button, adds to database
    def qaddDB(self):
        coption = self.DBsel.currentText()
        smimol = self.etDBsmi.text()
        sminame = self.etDBname.text()
        smident = self.DBdent.currentText()
        smicat = self.etDBsmicat.text()
        smigrps = self.etDBgrps.text()
        smictg = self.etDBctg.currentText()
        if self.lFFb.getState() and self.lFFa.getState():
            ffopt = 'BA'
        elif self.lFFb.getState() and not self.lFFa.getState():
            ffopt = 'B'
        elif not self.lFFb.getState() and self.lFFa.getState():
            ffopt = 'A'
        else:
            ffopt = 'N'
        if smictg == 'all':
            smictg = 'build functionalize'
        if smimol=='' or sminame=='':
            choice = QMessageBox.warning(self.DBWindow,'Error','Please specify molecule and name!')
        else:
            # add to database
            if 'ligand' in coption:
                emsg = addtoldb(smimol,sminame,smident,smicat,smigrps,smictg,ffopt)
            elif 'core' in coption:
                emsg = addtocdb(smimol,sminame,smicat)
            elif 'bind' in coption:
                emsg = addtobdb(smimol,sminame)
            if emsg:
                choice = QMessageBox.warning(self.DBWindow,'Error',emsg)
            else:
                choice = QMessageBox.information(self.DBWindow,'Add','Successfully added to the database!')
    ### callback for removal button, removes from db
    def qdelDB(self):
        coption = self.DBsel.currentIndex()
        sminame = self.etDBname.text()
        if sminame=='':
            choice = QMessageBox.warning(self.DBWindow,'Error','Please specify name!')
        else:
            # remove from database
            emsg = removefromDB(sminame,int(coption))
            if emsg:
                choice = QMessageBox.warning(self.DBWindow,'Error',emsg)
            else:
                choice = QMessageBox.information(self.DBWindow,'Remove','Removed from the database!')
    ### db type change
    def dbchange(self):
        ci = self.DBsel.currentIndex()
        if (ci==1):
            self.rtDBsmident.setDisabled(False)
            self.DBdent.setDisabled(False)
            self.rtDBsmicat.setDisabled(False)
            self.etDBsmicat.setDisabled(False)
            self.etDBgrps.setDisabled(False)
            self.etDBgrps.setDisabled(False)
            self.lFFb.setDisabled(False)
            self.lFFa.setDisabled(False)
        elif(ci==0):
            self.rtDBsmident.setDisabled(True)
            self.DBdent.setDisabled(True)
            self.rtDBsmicat.setDisabled(False)
            self.etDBsmicat.setDisabled(False)
            self.etDBgrps.setDisabled(True)
            self.etDBgrps.setDisabled(True)
            self.lFFb.setDisabled(True)
            self.lFFa.setDisabled(True)
        else:
            self.rtDBsmident.setDisabled(True)
            self.DBdent.setDisabled(True)
            self.rtDBsmicat.setDisabled(True)
            self.etDBsmicat.setDisabled(True)
            self.etDBgrps.setDisabled(True)
            self.etDBgrps.setDisabled(True)
            self.lFFb.setDisabled(True)
            self.lFFa.setDisabled(True)
    ### perform post-processing
    def postprocGUI(self):
        rdir = self.etpdir.text()
        if rdir[-1]=='/':
            rdir = rdir[:-1]
        defaultparams = ['main.py','-i',rdir+'/postproc.inp']
        grabguivarsP(self)
        self.pWindow.close()
        self.sgrid.setCurrentWidget(self.iWind)
        emsg = startgen(defaultparams,True,self)
        if not emsg:
            choice = QMessageBox.information(self.pWindow,'DONE','Your results are ready..')
    #############################
    #### Add geometry window ####
    #############################
    ### load molecule from file
    def qgeomload(self):
        name = QFileDialog.getOpenFileName(self.DBWindow,'Open File','.',"Molecule files *.xyz (*.xyz)")
        if name[0] != '':
            self.etgf.setText(os.path.relpath(name[0]))
    ### enable add new geometry
    def addgeom(self):
        self.geWindow.setWindowModality(2)
        self.geWindow.show()
        ### callback for addition button, adds to database
    def qaddg(self):
        globs = globalvars()
        gname = self.etgname.text().lower()
        gname = gname.replace(' ','_')
        gshort = self.etgshort.text()
        gfile = self.etgf.text()
        if gname=='' or gfile=='':
            choice = QMessageBox.warning(self.geWindow,'Error','Please specify name and xyz file!')
        elif not glob.glob(gfile):
            choice = QMessageBox.warning(self.geWindow,'Error','XYZ file '+gfile+' does not exist!')
        else:
            f = open(globs.installdir+'/Data/coordinations.dict','r')
            s = f.read().splitlines()
            f.close()
            if gname.lower() in s:
                choice = QMessageBox.warning(self.geWindow,'Add','Coordination '+gname+' already exists.')
                return
            # get geometry from xyz file
            f = open(gfile,'r')
            snew = f.read().splitlines()
            f.close()
            dent = int(snew[0])-1
            xyzl = ''
            for ii in range(0,dent+1):
                l = filter(None,re.split(' |\t',snew[2+ii]))
                xyzl += l[1]+' '+l[2]+' '+l[3]+'\n'
            # write new entry in coordinations.dict
            s.append(str(dent)+': '+gname+' '+gshort)
            ssort = filter(None,list(sorted(s[1:])))
            f = open(globs.installdir+'/Data/coordinations.dict','w')
            f.write(s[0]+'\n')
            for ss in ssort:
                f.write(ss+'\n')
            f.close()
            # write new backbone file
            f = open(globs.installdir+'/Data/'+gshort+'.dat','w')
            f.write(xyzl)
            f.close()
            choice = QMessageBox.information(self.geWindow,'Add','Successfully added to the database!')
        self.matchgeomcoord()
    ### callback for removal button, removes from db
    def qdelg(self):
        globs = globalvars()
        gname = self.etgname.text().lower()
        gname = gname.replace(' ','_')
        gshort = self.etgshort.text()
        if gname=='' and gshort=='':
            choice = QMessageBox.warning(self.geWindow,'Error','Please specify geometry name!')
        else:
            f = open(globs.installdir+'/Data/coordinations.dict','r')
            s = f.read()
            f.close()
            if gname.lower() not in s and gshort.lower() not in s:
                choice = QMessageBox.warning(self.geWindow,'Remove','Coordination '+gname+' does not exist.')
                return
            # remove entry from coordinations.dict
            snew = ''
            srem = ''
            for ss in s.splitlines():
                sl = filter(None,ss.split(' '))
                if gname.lower()!=sl[1] and gshort.lower()!=sl[2]:
                    snew += ss+'\n'
                else:
                    srem = filter(None,ss.split(' ')[-1])
            f = open(globs.installdir+'/Data/coordinations.dict','w')
            f.write(snew)
            f.close()
            # remove file
            if glob.glob(globs.installdir+'/Data/'+srem+'.dat'):
                os.remove(globs.installdir+'/Data/'+srem+'.dat')
            choice = QMessageBox.information(self.geWindow,'Remove','Successfully removed from the database!')
        self.matchgeomcoord()
    #############################
    ### Chem Database window ####
    #############################
    ### load molecule from file
    def qcDBload(self):
        name = QFileDialog.getOpenFileName(self.DBWindow,'Open File','.',"Molecule files *.xyz, *.mol *.sdf *.smi (*.xyz *.mol *.sdf *.smi)")
        if name[0] != '':
            self.etcDBsmi.setText(os.path.relpath(name[0]))
    ### enable add to database interface
    def searchDBW(self):
        globs = globalvars()
        self.cDBWindow.setWindowModality(2)
        self.cDBWindow.show()
        writef = False
        instdir = globs.installdir
        mwfn = globs.multiwfn
        cdbdir = globs.chemdbdir
        if not os.path.isdir(globs.chemdbdir):
            choice = QMessageBox.question(self.cDBWindow,'Database setup','It looks like the Chemical Database directory is not configured or does not exist. Would you like to configure it now?',
                QMessageBox.Yes, QMessageBox.No)
            if choice == QMessageBox.Yes:
                QMessageBox.information(self.cDBWindow,'Chem DB',"Please select the directory containing chemical databases.")
                cdbdir = QFileDialog.getExistingDirectory(self.cDBWindow,'Select the directory containing chemical databases.')
                if len(cdbdir) > 0:
                    writef = True
        if writef:
            f = open(globs.homedir+'/.molSimplify','w')
            f.write("INSTALLDIR="+instdir+'\n')
            f.write("CHEMDBDIR="+cdbdir+'\n')
            if len(mwfn) > 1 :
                f.write("MULTIWFN="+mwfn[0]+'\n')
            f.close()
            # get existing databases
            globsnew = globalvars()
            dbdir = globsnew.chemdbdir
            dbs0 = glob.glob(dbdir+"/*.sdf")
            dbs1 = [d.rsplit('/',1)[-1] for d in dbs0]
            dbs = [d.split('.',1)[0] for d in dbs1]
            for d in dbs:
                self.cDBsel.addItem(d)
    ### callback for database search
    def qaddcDB(self):
         ### collects all the info and passes it to molSimplify ###
        rdir = self.etrdir.text()
        if rdir[-1]=='/':
            rdir = rdir[:-1]
        # create running dir if not existing
        if not os.path.isdir(rdir):
            os.mkdir(rdir)
        args = grabdbguivars(self)
        defaultparams = ['main.py','-i',rdir+'/dbinput.inp']
        emsg = startgen(defaultparams,True,self)
        if not emsg:
            choice = QMessageBox.information(self.cDBWindow,'DONE','Search is done..')
        self.sgrid.setCurrentWidget(self.cDBWindow)
    ### db type change
    def cdbchange(self):
        ci = self.DBsel.currentIndex()
        if (ci==1):
            self.rtDBsmident.setDisabled(False)
            self.DBdent.setDisabled(False)
            self.rtDBsmicat.setDisabled(False)
            self.etDBsmicat.setDisabled(False)
        else:
            self.rtDBsmident.setDisabled(True)
            self.DBdent.setDisabled(True)
            self.rtDBsmicat.setDisabled(True)
            self.etDBsmicat.setDisabled(True)
    ####################
    ### Main window ####
    ####################
    ### run generation
    def runGUI(self):
        ### collects all the info and passes it to molSimplify ###
        rdir = self.etrdir.text()
        if rdir[-1]=='/':
            rdir = rdir[:-1]
        # create running dir if not existing
        if not os.path.isdir(rdir):
            try:
                os.mkdir(rdir)
            except:
                emsg = 'Directory '+rdir+' could not be created. Check your input.\n'
                QMessageBox.critical(self.wmain,'Problem',emsg)
                return
        # get parameters
        args = grabguivars(self)
        defaultparams = ['main.py','-i',rdir+'/geninput.inp']
        self.sgrid.setCurrentWidget(self.iWind)
        msgBox = QMessageBox()
        if self.randomchk.getState():
            msgBox.setText("Random generation initiated. This process might take some time.")
            msgBox.setIcon(1)
            msgBox.setInformativeText('Please be patient. OK?')
            msgBox.setWindowTitle('Running..')
            msgBox.exec_()
        # do the generation
        emsg = startgen(defaultparams,True,self)
        if not emsg:
            QMessageBox.information(self.wmain,'Done','Structure generation terminated successfully!')
        else:
            QMessageBox.warning(self.wmain,'Problem',emsg)
    ### draw ligands
    def drawligs(self):
        ### collects all the info and passes it to molSimplify ###
        rdir = self.etrdir.text()
        if rdir[-1]=='/':
            rdir = rdir[:-1]
        # creat running dir if not existing
        if not os.path.isdir(rdir):
            os.mkdir(rdir)
        args = grabguivars(self)
        if len(args['-lig']) < 1:
            qm = mQDialogWarn('Warning','No ligands are specified.')
            return False
        else:
            args['-lig']=args['-lig'].replace(' ','')
            lls = args['-lig'].split(',')
            liglist = []
            # check if multiple ligands in .smi file
            for l in lls:
                if '.smi' in l:
                    f = open(l,'r')
                    smis = filter(None,f.read().splitlines())
                    liglist += smis
                else:
                    liglist.append(l)
            globs = globalvars()
            licores = readdict(globs.installdir+'/Ligands/ligands.dict')
            ligs = []
            for l in liglist:
                if isinstance(l,unicode):
                    ll = unicodedata.normalize('NFKD',l).encode('ascii','ignore')
                else:
                    ll = l
                lig,emsg = lig_load(globs.installdir+'/',ll,licores)
                if emsg:
                    mQDialogWarn('Error',emsg)
                else:
                    ligs.append(lig.OBmol)
            if len(ligs)==0:
                return
            fcount = 0
            while glob.glob(rdir+'/ligs'+str(fcount)+'.png'):
                fcount += 1
            outputf  = 'ligs.smi'
            locf = 'ligs'+str(fcount)
            outbase = rdir+'/'+locf
            outf = pybel.Outputfile("smi",outputf,overwrite=True)
            for mol in ligs:
                outf.write(mol)
            # convert to svg
            if globs.osx:
                cmd = "/usr/local/bin/obabel -ismi "+outputf+" -O "+locf+".svg -xC -xi"
            else:
                cmd = "obabel -ismi "+outputf+" -O "+locf+".svg -xC -xi"
            t = mybash(cmd)
            print t
            if glob.glob(outputf):
                os.remove(outputf)
            else:
                mQDialogErr('Error','Image could not be generated\n.')
                return
            ####################
            ### draw ligands ###
            ####################
            if globs.osx:
                cmd = '/usr/local/bin/convert -density 1200 '+locf+'.svg '+locf+'.png'
            else:
                cmd = 'convert -density 1200 '+locf+'.svg '+locf+'.png'
            s = mybash(cmd)
            print s
            if not glob.glob(locf+'.png') :
                mQDialogInf('Done','2D representation of ligands generated in file ' +outbase+'.svg ! Conversion to png failed.')
            else:
                os.remove(locf+".svg")
                shutil.move(locf+'.png',outbase+'.png')
                self.c1p = mQPixmap(outbase+'.png')
                rows = self.lgrid.rowCount()
                if rows > 1:
                    for i in reversed(range(self.lgrid.count())): 
                        self.lgrid.itemAt(i).widget().setParent(None)
                self.lgrid.addWidget(self.c1p,0,0)
                # button for closing window
                ctip = 'Close current window'
                self.lwindow.setWindowTitle('Ligands 2D')
                self.lwclose = mQPushButton('Close',ctip,14)
                self.lwclose.clicked.connect(self.qcloseligs)
                self.lgrid.addWidget(self.lwclose,1,0)
                self.lwindow.showMaximized()
                center(self.lwindow)
    ### draw ligands
    def viewgeom(self):
        globs = globalvars()
        # get geometry
        geom = self.dcoordg.currentText()
        gfname = globs.installdir+'/icons/geoms/'+geom+'.png'
        if glob.glob(gfname):
            rows = self.lgrid.rowCount()
            if rows > 1:
                for i in reversed(range(self.lgrid.count())): 
                    self.lgrid.itemAt(i).widget().setParent(None)
            self.c1p = mQPixmap(gfname)
            self.lgrid.addWidget(self.c1p,0,0)
            # button for closing window
            ctip = 'Close current window'
            self.lwclose = mQPushButton('Close',ctip,14)
            self.lwclose.clicked.connect(self.qcloseligs)
            self.lgrid.addWidget(self.lwclose,1,0)
            self.lwindow.setWindowTitle('Geometry:'+self.dcoordg.currentText())
            self.lwindow.show()
            center(self.lwindow)
        else:
            mQDialogWarn('Warning','No file '+gfname+' exists..')
    ### draw results from db search
    def drawres(self):
        ### collects all the info and passes it to molSimplify ###
        rdir = self.etrdir.text()
        if rdir[-1]=='/':
            rdir = rdir[:-1]
        # creat running dir if not existing
        if not os.path.isdir(rdir):
            os.mkdir(rdir)
        outf = rdir+'/'+self.etcDBoutf.text()
        outf = outf.replace(' ','')+'.smi'
        if not glob.glob(outf):
            mQDialogWarn('Warning','No database results file in '+rdir)
            return False
        else:
            lls = [outf]
            liglist = []
            # check if multiple ligands in .smi file
            for l in lls:
                if '.smi' in l:
                    f = open(l,'r')
                    smis = filter(None,f.read().splitlines())
                    liglist += smis
                else:
                    liglist.append(l)
            globs = globalvars()
            licores = readdict(globs.installdir+'/Ligands/ligands.dict')
            ligs = []
            for l in liglist:
                if isinstance(l,unicode):
                    ll = unicodedata.normalize('NFKD',l).encode('ascii','ignore')
                else:
                    ll = l
                lig,emsg = lig_load(globs.installdir+'/',ll,licores)
                if not emsg:
                    ligs.append(lig.OBmol)
            if len(ligs)==0:
                return
            fcount = 0
            while glob.glob(rdir+'/ligs'+str(fcount)+'.png'):
                fcount += 1
            outputf  = 'ligs.smi'
            locf = 'ligs'+str(fcount)
            outbase = rdir+'/'+locf
            outf = pybel.Outputfile("smi",outputf,overwrite=True)
            for mol in ligs:
                outf.write(mol)
            # convert to svg
            if globs.osx:
                cmd = "/usr/local/bin/obabel -ismi "+outputf+" -O "+locf+".svg -xC -xi"
            else:
                cmd = "obabel -ismi "+outputf+" -O "+locf+".svg -xC -xi"
            t = mybash(cmd)
            print t
            if glob.glob(outputf):
                os.remove(outputf)
            else:
                mQDialogInf('Error','Image could not be generated\n.')
                return
            ####################
            ### draw ligands ###
            ####################
            if globs.osx:
                cmd = '/usr/local/bin/convert -density 1200 '+locf+'.svg '+locf+'.png'
            else:
                cmd = 'convert -density 1200 '+locf+'.svg '+locf+'.png'
            s = mybash(cmd)
            print s
            if not glob.glob(locf+'.png') :
                mQDialogInf('Done','2D representation of ligands generated in file ' +outbase+'.svg ! Conversion to png failed.')
            else:
                os.remove(locf+".svg")
                shutil.move(locf+'.png',outbase+'.png')
                self.c1p = mQPixmap(outbase+'.png')
                rows = self.lgrid.rowCount()
                if rows > 1:
                    for i in reversed(range(self.lgrid.count())): 
                        self.lgrid.itemAt(i).widget().setParent(None)
                self.lgrid.addWidget(self.c1p,0,0)
                # button for closing window
                ctip = 'Close current window'
                self.lwclose = mQPushButton('Close',ctip,14)
                self.lwclose.clicked.connect(self.qcloseligs)
                self.lgrid.addWidget(self.lwclose,1,0)
                self.lwindow.showMaximized()
                center(self.lwindow)
    ### enable random input
    def enablerandom(self):
        if self.randomchk.isChecked():
            self.rtrgen.setDisabled(False)
            self.etrgen.setDisabled(False)
            self.rtlignum.setDisabled(False)
            self.etlignum.setDisabled(False)
            self.rtliggrp.setDisabled(False)
            self.etliggrp.setDisabled(False)
            self.rtligctg.setDisabled(False)
            self.etligctg.setDisabled(False)
            self.randkHs.setDisabled(False)
            if len(self.etrgen.text())==0:
                self.etrgen.setText('1')
        else:
            self.rtrgen.setDisabled(True)
            self.etrgen.setDisabled(True)
            self.rtlignum.setDisabled(True)
            self.etlignum.setDisabled(True)
            self.rtliggrp.setDisabled(True)
            self.etliggrp.setDisabled(True)
            self.rtligctg.setDisabled(True)
            self.etligctg.setDisabled(True)
            self.randkHs.setDisabled(True)
    ### generate all enable FF input
    def disableffinput(self):
        if self.chkgenall.isChecked():
            self.dff.setDisabled(True)
            self.dffba.setDisabled(True)
            self.chkFF.state = True
            self.chkFF.setDisabled(True)
        else:
            self.dff.setDisabled(False)
            self.dffba.setDisabled(False)
            self.chkFF.setDisabled(False)
            self.enableffinput()
    ### enable FF input
    def enableffinput(self):
        if self.chkFF.isChecked():
            self.chkFF.state = False
            self.dff.setDisabled(False)
            self.dffba.setDisabled(False)
        else:
            self.dff.setDisabled(True)
            self.dffba.setDisabled(True)
            self.chkFF.state = True
    ### load file for editing
    def qloadinput(self):
        name = QFileDialog.getOpenFileName(self.wmain,'Open File')[0]
        if name != '':
            loadfrominputfile(self,name)
    ### load directory
    def dirload(self):
        name = QFileDialog.getExistingDirectory(self.wmain,'Select Directory')
        if len(name) >0 and name[0] != '':
            self.etrdir.setText(name)
    ### save as input file
    def qsaveinput(self,gui):
        name = QFileDialog.getSaveFileName(self.wmain,'Save as..','.',"Input files * (*)")[0]
        if name != '':
            varsg = grabguivars(self)
            writeinputc(varsg,name)
    ### show help menu
    def qshowhelp(self):
        globs = globalvars()
        QMessageBox.information(self.wmain,'About',globs.about)
    def getscreensize(self):
        screenShape = QDesktopWidget().screenGeometry()
        width = int(screenShape.width())
        height = int(screenShape.height())
        return [width,height]
    ### slider changed value
    def sliderChanged(self,val):
        self.distper.setText('Distort:'+str(val)+'%')
    ### match index with coordination
    def matchgeomcoord(self):
        globs = globalvars()
        # get current index
        dc=self.dcoord.currentIndex()
        coords,geomnames,geomshorts,geomgroups = getgeoms()
        qcav = geomgroups
        ctip = geomnames
        # empty the box
        for i in range(0,self.dcoordg.count()):
            self.dcoordg.removeItem(0)
        qc = qcav[dc]
        # add to box
        for i,t in enumerate(qc):
            self.dcoordg.addItem(QIcon(globs.installdir+'/icons/geoms/'+t+'.png'),t)
        self.dcoordg.setIconSize(QSize(60,60))
        # set default geometry
        self.dcoordg.setCurrentIndex(0)
        # get global index
        elem = [i for i,s in enumerate(geomshorts) if s in self.dcoordg.currentText()]
        ct = ''
        for ii in range(elem[0],elem[0]+len(geomgroups[dc])):
            ct += ctip[ii]+', '
        ct = ct[:-2]
        # set correct tooltip
        self.dcoordg.setToolTip(ct)
    ### enable extra molecule input
    def enableemol(self):
        if self.chkM.isChecked():
            self.chkM.state = False
            self.txtamol.setDisabled(False)
            self.rtbind.setDisabled(False)
            self.etbind.setDisabled(False)
            self.rtbsmi.setDisabled(False)
            self.etbsmi.setDisabled(False)
            self.rtnbind.setDisabled(False)
            self.etnbind.setDisabled(False)    
            self.rtplace.setDisabled(False)
            self.etplacemin.setDisabled(False)
            self.etplacemax.setDisabled(False)
            self.dmolp.setDisabled(False)
            self.rtchbind.setDisabled(False)
            self.etchbind.setDisabled(False)
            self.rtplacea.setDisabled(False)
            self.etplacephi.setDisabled(False)
            self.etplacetheta.setDisabled(False)
            #self.rtmaskbind.setDisabled(False)
            self.etmaskbind.setDisabled(False)
            self.chsep.setDisabled(False)
        else:
            self.txtamol.setDisabled(True)
            self.rtbind.setDisabled(True)
            self.etbind.setDisabled(True)
            self.rtbsmi.setDisabled(True)
            self.etbsmi.setDisabled(True)
            self.rtnbind.setDisabled(True)
            self.etnbind.setDisabled(True)    
            self.rtplace.setDisabled(True)
            self.etplacemin.setDisabled(True)
            self.etplacemax.setDisabled(True)
            self.dmolp.setDisabled(True)
            self.chkM.state = True
            self.rtchbind.setDisabled(True)
            self.etchbind.setDisabled(True)
            self.rtplacea.setDisabled(True)
            self.etplacephi.setDisabled(True)
            self.etplacetheta.setDisabled(True)
            #self.rtmaskbind.setDisabled(True)
            self.etmaskbind.setDisabled(True)
            self.chsep.setDisabled(True)
    #####################
    #### QEt/g input ####
    #####################
    ### callback for QE input
    def qcinput(self):
        if self.chch.getState():
            self.etqcgch.setDisabled(True)
            self.etqctch.setDisabled(True)
            self.etqcQch.setDisabled(True)
        if self.qcode.currentIndex()==0: # generate terachem input
            self.qctWindow.setWindowModality(2)
            self.qctWindow.show()
        elif self.qcode.currentIndex()==1: # generate GAMESS input
            self.qcgWindow.setWindowModality(2)
            self.qcgWindow.show()
        elif self.qcode.currentIndex()==2: # generate Qchem input
            self.qcQWindow.setWindowModality(2)
            self.qcQWindow.show()
    # make default button callback
    def jobdef(self):
        grabguivarsjob(self)
        QMessageBox.information(self.wmain,'Done','The current settings are the default ones now.')
    def qctdef(self):
        grabguivarstc(self)
        QMessageBox.information(self.wmain,'Done','The current settings are the default ones now.')
    def qcgdef(self):
        grabguivarsgam(self)
        QMessageBox.information(self.wmain,'Done','The current settings are the default ones now.')
    def qcqdef(self):
        grabguivarsqch(self)
        QMessageBox.information(self.wmain,'Done','The current settings are the default ones now.')
    def qcgload(self):
        name = QFileDialog.getOpenFileName(self.qctWindow,'Open File','.',"GAMESS input files")
        if name[0] != '':
            f = open(name[0],'r')
            self.qcgWindow.molf = f.read()
            f.close()
    ### enable QE input
    def enableqeinput(self):
        if self.chkI.isChecked():
            self.chkI.state = False
            self.qcode.setDisabled(False)
            self.butQc.setDisabled(False)
            self.chch.setDisabled(False)
        else:
            self.qcode.setDisabled(True)
            self.butQc.setDisabled(True)
            self.chch.setDisabled(True)
            self.chkI.state = True
    #########################
    #### jobscript input ####
    #########################
    ### enable Jobscript input
    def enablejinput(self):
        if self.chkJ.isChecked():
            self.chkJ.state = False
            self.scheduler.setDisabled(False)
            self.butJob.setDisabled(False)
        else:
            self.scheduler.setDisabled(True)
            self.butJob.setDisabled(True)
            self.chkJ.state = True
    def jobenable(self):
            self.jWindow.setWindowModality(2)
            self.jWindow.show()
    ### load file
    def jload(self):
        name = QFileDialog.getOpenFileName(self.qctWindow,'Open File','.',"Jobscript files")
        if name[0] != '':
            f = open(name[0],'r')
            self.jWindow.molf = f.read()
            f.close()
    ###############################
    #### post-processing input ####
    ###############################
    ### enable Post processing input
    def setupp(self):
        self.pWindow.setWindowModality(2)
        self.pWindow.show()
        # check if Multiwfn exists
        globs = globalvars()
        inputtxt = '0\n0\n' 
        f = open('input1','w')
        f.write(inputtxt)
        f.close()
        writef = False
        instdir = globs.installdir
        mwfn = globs.multiwfn
        cdbdir = globs.chemdbdir
        if not os.path.isfile(globs.multiwfn[1:-1]):
            choice = QMessageBox.question(self.pWindow,'Multiwfn setup','It looks like the Multiwfn executable is not configured or does not exist. Would you like to configure it now?',
                QMessageBox.Yes, QMessageBox.No)
            if choice == QMessageBox.Yes:
                QMessageBox.information(self.pWindow,'Multiwfn',"Please select the Multiwfn executable.")
                mwfn = QFileDialog.getOpenFileName(self.pWindow,'Select the Multiwfn executable.','.')
                if len(mwfn[0]) > 1:
                    writef = True
        if writef:
            f = open(globs.homedir+'/.molSimplify','w')
            f.write("INSTALLDIR="+instdir+'\n')
            if len(cdbdir) > 0:
                f.write("CHEMDBDIR="+cdbdir+'\n')
            if len(mwfn) > 0 :
                f.write("MULTIWFN="+mwfn[0]+"\n")
            f.close()
        newglobs = globalvars()
        com = newglobs.multiwfn
        tt = mybash(com + '< input1')
        os.remove('input1')
        if not 'Multifunctional Wavefunction Analyzer' in tt:
            self.pch.setDisabled(True)
            self.pwfnav.setDisabled(True)
            self.pcub.setDisabled(True)
            self.pdeloc.setDisabled(True)
    ### load directory
    def pdload(self):
        name = QFileDialog.getExistingDirectory(self.pWindow,'Select Directory')
        if len(name) >0 and name[0] != '':
            self.etpdir.setText(name)
    ###############################
    ###### general callbacks ######
    ###############################
    ### exit application ###
    def qexit(self):
        choice = QMessageBox.question(self.wmain,'Exit','Are you sure you want to quit?',
                QMessageBox.Yes, QMessageBox.No)
        if choice == QMessageBox.Yes:
            sys.exit()
        else:
            pass
    ### hide current widget ###
    def qhide(self):
        self.hide()
    ### close ligands window
    def qcloseligs(self):
        self.lwindow.hide()
    ### return to main window ###
    def qretmain(self):
        # hide all windows
        self.qctWindow.hide()
        self.qcgWindow.hide()
        self.qcQWindow.hide()
        self.pWindow.hide()
        self.jWindow.hide()
        self.DBWindow.hide()
        self.cDBWindow.hide()
        self.geWindow.hide()

