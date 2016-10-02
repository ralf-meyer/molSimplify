# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

##############################################################
########## This script processes the input file  #############
##############################################################

# import std modules
import glob, os, re, argparse, sys
from io import *
from Classes.globalvars import *

######################################################
########## check core/ligands specified  #############
######################################################
### checks input for correctness ###
def checkinput(args):
    emsg = False
    # check if core is specified
    if not (args.core):
            emsg = 'You need to specify at least the core of the structures.\n'
            print emsg
            return emsg
    # check if ligands are specified
    if not args.lig and not args.rgen:
        if args.gui:
            from Classes.mWidgets import mQDialogWarn
            qqb = mQDialogWarn('Warning','You specified no ligands. Please use the -lig flag. Core only generation..')
            qqb.setParent(args.gui.wmain)
        else:
            print 'WARNING: You specified no ligands. Please use the -lig flag. Forced generation..\n'
    return emsg

###########################################
########## check true or false  ###########
###########################################
def checkTrue(arg):
    if 'y' in arg.lower() or '1' in arg.lower() or 't' in arg.lower() or arg==1:
        return True
    else:
        return False
                
            
###########################################
########## consolidate lists  #############
###########################################
### consolidate arguments ###
def cleaninput(args):
    globs = globalvars()
    # check ligands
    if args.lig:
        ls = []
        ligdic = readdict(globs.installdir+'/Ligands/simple_ligands.dict')
        for i,s in enumerate(args.lig):
            if isinstance(s,list):
                for ss in s:
                    if ss in ligdic.keys():
                        ss = ligdic[ss][0]
                    ls.append(ss)
            else:
                if s in ligdic.keys():
                        s = ligdic[s][0]
                ls.append(s)
        args.lig = ls
    # check sminame
    if args.sminame:
        ls = []
        for i,s in enumerate(args.sminame):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.sminame = ls
    # check qoption
    if args.qoption:
        ls = []
        for i,s in enumerate(args.qoption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.qoption = ls
    # check sysoption
    if args.sysoption:
        ls = []
        for i,s in enumerate(args.sysoption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.sysoption = ls
    # check ctrloption
    if args.ctrloption:
        ls = []
        for i,s in enumerate(args.ctrloption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.ctrloption = ls
    # check scfoption
    if args.scfoption:
        ls = []
        for i,s in enumerate(args.scfoption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.scfoption = ls
    # check statoption
    if args.statoption:
        ls = []
        for i,s in enumerate(args.statoption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.statoption = ls
    # check remoption
    if args.remoption:
        ls = []
        for i,s in enumerate(args.remoption):
            if isinstance(s,list):
                for ss in s:
                    ls.append(ss)
            else:
                ls.append(s)
        args.remoption = ls
    # convert keepHs to boolean
    if args.keepHs:
        for i,s in enumerate(args.keepHs):
            args.keepHs[i]=checkTrue(s)
    # convert ff option to abe code
    if args.ff and args.ffoption:
        b = False
        a = False
        e = False
        opts = args.ffoption
        args.ffoption = ''
        for op in opts:
            op = op.strip(' ')
            if op[0].lower()=='b':
                args.ffoption += 'b'
            if op[0].lower()=='a':
                args.ffoption += 'a'
    elif args.ff:
        args.ffoption = 'ba'

###################################################
##########  parse command line input  #############
###################################################
### parses inputfile ###
def parseCLI(args):
    cliargs = ' '.join(args)
    s = filter(None,cliargs.split('-'))
    fname = 'CLIinput.inp'
    f = open(fname,'w')
    f.write('# molSimplify input file generated from CLI input\n')
    for line in s:
       f.write('-'+line+'\n')
    f.close()
    return fname
    
###########################################
##########  parse input file  #############
###########################################
### parses inputfile ###
def parseinput(args):
    for line in open(args.i):
        if '-lig' not in line and '-core' not in line and '-bind' not in line:
            line = line.split('#')[0] # remove comments
        li = line.strip()
        if not li.startswith("#") and len(li)>0: # remove comments/empty lines
            l = li.split('\n')[0] # remove newlines
            l = filter(None,re.split(' |,|\t|&',l))
            # parse general arguments
            if (l[0]=='-core'):
                args.core = [ll for ll in l[1:]]
            if (l[0]=='-ccatoms'):
                args.ccatoms = [int(ll)-1 for ll in l[1:]]
            if (l[0]=='-rundir'):
                args.rundir = line.split("#")[0].strip('\n')
                args.rundir = args.rundir.split('-rundir')[1]
                args.rundir = args.rundir.lstrip(' ')
                if (args.rundir[-1]=='/'):
                    args.rundir = args.rundir[:-1]
            if (l[0]=='-suff'):
                args.suff = l[1].strip('\n')
            ### parse structure generation arguments ###
            if (l[0]=='-bind'):
                l = filter(None,re.split(' |,|\t',line[:-1]))
                # discard comments
                for ibind,lbind in enumerate(l):
                    if lbind=='#':
                        l=l[:ibind]
                        break
                args.bind = l[1]
            if (l[0]=='-nbind'):
                args.bindnum = l[1]
            if (l[0]=='-bcharge'):  # parse charge for binding species
                args.bcharge = l[1]
            if (l[0]=='-btheta'):
                args.btheta = l[1]
            if (l[0]=='-bphi'):
                args.bphi = l[1]
            if (l[0]=='-bsep'):
                args.bsep = l[1]
            if (l[0]=='-bref'):
                args.bref = l[1:]
            if (l[0]=='-nambsmi'):
                args.nambsmi = l[1]
            if (l[0]=='-maxd'):
                args.maxd = l[1]
            if (l[0]=='-mind'):
                args.mind = l[1]
            if (l[0]=='-oxstate'):
                args.oxstate = l[1]
            if (l[0]=='-coord'):
                args.coord = l[1]
            if (l[0]=='-geometry'):
                args.geometry = l[1].lower()
            # parse ligands
            if (l[0]=='-lig'):
                l = filter(None,re.split(' |,|\t',line[:-1]))
                # discard comments
                for ilig,llig in enumerate(l):
                    if llig=='#':
                        l=l[:ilig]
                        break
                if args.lig:
                    args.lig.append(l[1:])
                else:
                    args.lig = l[1:]
            if (l[0]=='-lignum'):
                args.lignum = l[1]
            if (l[0]=='-liggrp'):
                args.liggrp = l[1]
            if (l[0]=='-ligctg'):
                args.ligctg = l[1]
            if (l[0]=='-ligocc'):
                args.ligocc = l[1:]
            if (l[0]=='-rkHs'):
                args.rkHs = checkTrue(l[1])
            if (l[0]=='-ligloc'):
                args.ligloc = checkTrue(l[1])
            if (l[0]=='-ligalign'):
                args.ligalign = checkTrue(l[1])
            if (l[0]=='-replig'):
                args.replig = checkTrue(l[1])
            if (l[0]=='-genall'):
                args.genall = checkTrue(l[1])
            if (l[0]=='-MLbonds'):
                args.MLbonds = l[1:]
            if (l[0]=='-distort'):
                args.distort = l[1]
            if (l[0]=='-langles'):
                args.langles = l[1:]
            if (l[0]=='-pangles'):
                args.pangles = l[1:]
            if (l[0]=='-rgen'):
                args.rgen = l[1:]
            if (l[0]=='-keepHs'):
                args.keepHs = l[1:]
            if (l[0]=='-ff'):
                args.ff = l[1].lower()
            if (l[0]=='-ffoption'):
                args.ffoption = l[1:]
            if (l[0]=='-place'):
                args.place = l[1]
            if (l[0]=='-sminame'):
                if args.sminame:
                    args.sminame.append(l[1:])
                else:
                    args.sminame = l[1:]
            if 'smicat' in line:
                args.smicat = []
                l = line.split(' ',1)[1]
                l = filter(None,re.split('/|\t|&|\n',l))
                for ll in l:
                    lloc = []
                    l1 = filter(None,re.split(',| ',ll))
                    for lll in l1:
                        lloc.append(int(lll)-1)
                    args.smicat.append(lloc)
            # parse qc arguments
            if (l[0]=='-qccode'):
                args.qccode = l[1]
            if (l[0]=='-calccharge'):
                args.calccharge = checkTrue(l[1])
            if (l[0]=='-charge'):
                args.charge = l[1:]
            if (l[0]=='-spin'):
                args.spin = l[1:]
            if (l[0]=='-runtyp'):
                args.runtyp = l[1]
            if (l[0]=='-method'):
                args.method = l[1:]
            # parse terachem arguments
            if (l[0]=='-basis'):
                args.basis = l[1]
            if (l[0]=='-dispersion'):
                args.dispersion = l[1].strip('\n').lower()
            if (l[0]=='-qoption'):
                if args.qoption:
                    args.qoption.append(l[1:])
                else:
                    args.qoption = l[1:]
            # parse qchem arguments
            if (l[0]=='-exchange'):
                args.exchange = l[1]
            if (l[0]=='-correlation'):
                args.correlation = l[1]
            if (l[0]=='-unrestricted'):
                args.unrestricted = checkTrue(l[1])
            if (l[0]=='-remoption'):
                if args.remoption:
                    args.remoption.append(l[1:])
                else:
                    args.remoption = l[1:]
            # parse gamess arguments
            if (l[0]=='-gbasis'):
                args.gbasis = l[1]
            if (l[0]=='-ngauss'):
                args.ngauss = l[1]
            if (l[0]=='-npfunc'):
                args.npfunc = l[1]
            if (l[0]=='-ndfunc'):
                args.ndfunc = l[1]
            if (l[0]=='-sysoption'):
                if args.sysoption:
                    args.sysoption.append(l[1:])
                else:
                    args.sysoption = l[1:]
            if (l[0]=='-ctrloption'):
                if args.ctrloption:
                    args.ctrloption.append(l[1:])
                else:
                    args.ctrloption = l[1:]
            if (l[0]=='-scfoption'):
                if args.scfoption:
                    args.scfoption.append(l[1:])
                else:
                    args.scfoption = l[1:]
            if (l[0]=='-statoption'):
                if args.statoption:
                    args.statoption.append(l[1:])
                else:
                    args.statoption = l[1:]
            # parse jobscript arguments
            if (l[0]=='-jsched'):
                args.jsched = l[1]
            if (l[0]=='-jname'):
                args.jname = l[1]
            if (l[0]=='-memory'):
                args.memory = l[1]
            if (l[0]=='-wtime'):
                args.wtime = l[1]
            if (l[0]=='-queue'):
                args.queue = l[1]
            if (l[0]=='-gpus'):
                args.gpus = l[1]
            if (l[0]=='-cpus'):
                args.cpus = l[1]
            if (l[0]=='-modules'):
                args.modules = l[1:]
            if (l[0]=='-joption'):
                if not args.joption:
                    args.joption = []
                opts = ''
                for ll in l[1:]:
                    opts += ll+' '
                args.joption.append(opts)
            if (l[0]=='-jcommand'):
                if not args.jcommand:
                    args.jcommand = []
                opts = ''
                for ll in l[1:]:
                    opts += ll+' '
                args.jcommand.append(opts)
            # parse database arguments
            if (l[0]=='-dbsim'):
                args.dbsearch = True
                args.dbsim = l[1]
            if (l[0]=='-dbresults'):
                args.dbresults = l[1]
            if (l[0]=='-dboutputf'):
                args.dboutputf = l[1]
            if (l[0]=='-dbbase'):
                args.dbbase = l[1]
            if (l[0]=='-dbsmarts'):
                args.dbsearch = True
                args.dbsmarts = l[1]
            if (l[0]=='-dbcatoms'):
                args.dbcatoms = l[1:]
            if (l[0]=='-dbfinger'):
                args.dbfinger = l[1]
            if (l[0]=='-dbatoms'):
                args.dbatoms = l[1]
            if (l[0]=='-dbbonds'):
                args.dbbonds = l[1]
            if (l[0]=='-dbarbonds'):
                args.dbarbonds = l[1]
            if (l[0]=='-dbsbonds'):
                args.dbsbonds = l[1]
            if (l[0]=='-dbmw'):
                args.dbmw = l[1]
            # parse postprocessing arguments
            if (l[0]=='-postp'):
                args.postp = True
            if (l[0]=='-postqc'):
                args.postqc = l[1]
            if (l[0]=='-postdir'):
                args.postdir = line.split("#")[0].strip('\n')
                args.postdir = args.postdir.split('-postdir')[1]
                args.postdir = args.postdir.lstrip(' ')
            if (l[0]=='-pres'):
                args.pres = True
            if (l[0]=='-pwfninfo'):
                args.pwfninfo = True
            if (l[0]=='-pcharge'):
                args.pcharge = True
            if (l[0]=='-pgencubes'):
                args.pgencubes = True
            if (l[0]=='-porbinfo'):
                args.porbinfo = True
            if (l[0]=='-pdeloc'):
                args.pdeloc = True
            #if (l[0]=='-pdorbs'):
            #    args.pdorbs = True
            if (l[0]=='-pnbo'):
                args.pnbo = True
                
#############################################################
########## mainly for help and listing options  #############
#############################################################
### parses commandline arguments and prints help information ###
def parsecommandline(parser):
    globs = globalvars()
    installdir = globs.installdir+'/'
    # first variable is the flag, second is the variable in the structure. e.g -i, --infile assigns something to args.infile
    parser.add_argument("-i","--i",help="input file")
    # top directory options
    parser.add_argument("-rundir","--rundir",help="directory for jobs",action="store_true")
    parser.add_argument("-suff","--suff", help="suffix for jobs folder.",action="store_true")
    # structure generation options
    parser.add_argument("-ccatoms","--ccatoms", help="core connection atoms indices, indexing starting from 1",action="store_true")
    parser.add_argument("-coord","--coord", help="coordination such as 4,5,6",action="store_true") # coordination e.g. 6 
    parser.add_argument("-core","--core", help="core structure with currently available: "+getcores(installdir),action="store_true") #e.g. ferrocene
    parser.add_argument("-bind","--bind", help="binding species with currently available: "+getbinds(installdir),action="store_true") #e.g. bisulfate, nitrate, perchlorate -> For binding
    parser.add_argument("-bcharge","--bcharge",default='0', help="binding species charge, default 0",action="store_true") 
    parser.add_argument("-bphi","--bphi", help="azimuthal angle phi for binding species, default random between 0 and 180",action="store_true") 
    parser.add_argument("-bref","--bref", help="reference atoms for placement of extra molecules, default COM (center of mass). e.g. 1,5 or 1-5, Fe, COM",action="store_true") 
    parser.add_argument("-bsep","--bsep", help="flag for separating extra molecule in input or xyz file",action="store_true")
    parser.add_argument("-btheta","--btheta", help="polar angle theta for binding species, default random between 0 and 360",action="store_true") 
    parser.add_argument("-geometry","--geometry", help="geometry such as TBP (trigonal bipyramidal)",action="store_true") # geometry
    parser.add_argument("-genall","--genall", help="Generate complex both with and without FF opt.",action="store_true") # geometry
    parser.add_argument("-lig","--lig", help="ligand structure name or SMILES with currently available: "+getligs(installdir),action="store_true") #e.g. acetate (in smilesdict)
    parser.add_argument("-ligocc","--ligocc", help="number of corresponding ligands e.g. 2,2,1",action="store_true") # e.g. 1,2,1
    parser.add_argument("-lignum","--lignum", help="number of ligand types e.g. 2",action="store_true") 
    parser.add_argument("-liggrp","--liggrp", help="ligand group for random generation",action="store_true") 
    parser.add_argument("-ligctg","--ligctg", help="ligand category for random generation",action="store_true") 
    parser.add_argument("-rkHs","--rkHs", help="keep Hydrogens for random generation",action="store_true") 
    parser.add_argument("-ligloc","--ligloc", help="force location of ligands in the structure generation yes/True/1 or no/False/0",action="store_true") 
    parser.add_argument("-ligalign","--ligalign", help="smart alignment of ligands in the structure generation yes/True/1 or no/False/0",action="store_true") 
    parser.add_argument("-MLbonds","--MLbonds", help="custom M-L bond length for corresponding ligand in A e.g. 1.4",action="store_true") 
    parser.add_argument("-distort","--distort", help="randomly distort backbone. Ranges from 0 (no distortion) to 100. e.g. 20",action="store_true") 
    parser.add_argument("-langles","--langles", help="custom angles (polar theta, azimuthal phi) for corresponding ligand in degrees separated by '/' e.g. 20/30,10/20",action="store_true") 
    parser.add_argument("-pangles","--pangles", help="custom angles (polar theta, azimuthal phi) for corresponding connectino points in degrees separated by '/' e.g. 20/30,10/20",action="store_true") 
    parser.add_argument("-nbind","--bindnum", help="number of binding species copies for random placement",action="store_true") #different geometric arrangements for calculating binding energy
    parser.add_argument("-rgen","--rgen", help="number of random generated molecules",action="store_true")
    parser.add_argument("-replig","--replig", help="flag for replacing ligand at specified connection point",action="store_true")
    parser.add_argument("-ff","--ff",help="select force field for FF optimization. Available: MMFF94, UFF, GAFF, Ghemical",action="store_true")
    parser.add_argument("-ffoption","--ffoption",help="select when to perform FF optimization. Options: B(Before),A(After),BA",action="store_true")
    parser.add_argument("-keepHs","--keepHs", help="force keep hydrogens. By default ligands are stripped one hydrogen in order to connect to the core",action="store_true") 
    parser.add_argument("-smicat","--smicat", help="connecting atoms corresponding to smiles. Indexing starts at 1 which is the default value as well",action="store_true")
    parser.add_argument("-sminame","--sminame", help="name for smiles species used in the folder naming. e.g. amm",action="store_true") 
    parser.add_argument("-nambsmi","--nambsmi", help="name of SMILES string for binding species e.g. carbonmonoxide",action="store_true")
    parser.add_argument("-maxd","--maxd", help="maximum distance above cluster size for molecules placement maxdist=size1+size2+maxd", action="store_true")
    parser.add_argument("-mind","--mind", help="minimum distance above cluster size for molecules placement mindist=size1+size2+mind", action="store_true")
    parser.add_argument("-place","--place", help="place binding species relative to core. Takes either angle (0-360) or ax/s for axial side",action="store_true")
    parser.add_argument("-oxstate","--oxstate", help="oxidation state of the metal, used for bond lengths",action="store_true")
    # quantum chemistry options
    parser.add_argument("-qccode","--qccode", help="quantum chemistry code. Choices: TeraChem or GAMESS or QChem",action="store_true") 
    parser.add_argument("-charge","--charge", help="charge for system (default: neutral).",action="store_true")
    parser.add_argument("-calccharge","--calccharge", help="Flag to calculate charge.",action="store_true")
    parser.add_argument("-spin","--spin", help="spin multiplicity for system (default: singlet) e.g. 1",action="store_true")
    parser.add_argument("-runtyp","--runtyp", help="run type. Choices: optimization, energy",action="store_true")
    parser.add_argument("-method","--method", help="electronic structure method. Specify UDFT for unrestricted calculation(default: b3lyp) e.g. ub3lyp",action="store_true")
    # terachem arguments
    parser.add_argument("-basis","--basis", help="basis for terachem or qchem job (default: LACVP* or lanl2dz)",action="store_true")
    parser.add_argument("-dispersion","--dispersion", help="dispersion forces. Default: no e.g. d2,d3",action="store_true")
    parser.add_argument("-qoption","--qoption", help="extra arguments for TeraChem in syntax keyword value, e.g. maxit 100",action="store_true")
    # qchem arguments
    parser.add_argument("-exchange","--exchange",help="exchange in qchem job (default b3lyp)",action="store_true")
    parser.add_argument("-correlation","--correlation",help="correlation in qchem job (default none)",action="store_true")
    parser.add_argument("-remoption","--remoption", help="extra arguments for qchem $rem block in syntax keyword value, e.g. INCFOCK 0",action="store_true")
    parser.add_argument("-unrestricted","--unrestricted", help="unrestricted calculation, values: 0/1 False/True",action="store_true")
    # gamess arguments
    parser.add_argument("-gbasis","--gbasis", help="GBASIS option in GAMESS e.g. CCT",action="store_true")
    parser.add_argument("-ngauss","--ngauss", help="NGAUSS option in GAMESS e.g. N31",action="store_true")
    parser.add_argument("-npfunc","--npfunc", help="NPFUNC option for diffuse functions in GAMESS e.g. 2",action="store_true")
    parser.add_argument("-ndfunc","--ndfunc", help="NDFUNC option for diffuse functions in GAMESS e.g. 1",action="store_true")
    parser.add_argument("-sysoption","--sysoption", help="extra arguments for $SYSTEM GAMESS block in syntax keyword value, e.g. MWORDS 20",action="store_true")
    parser.add_argument("-ctrloption","--ctrloption", help="extra arguments for $CONTRL GAMESS block in syntax keyword value, e.g. ISPHER 1",action="store_true")
    parser.add_argument("-scfoption","--scfoption", help="extra arguments for $SCF GAMESS block in syntax keyword value, e.g. DIIS .TRUE.",action="store_true")
    parser.add_argument("-statoption","--statoption", help="extra arguments for $STATPT GAMESS block in syntax keyword value, e.g. NSTEP 100",action="store_true")
    # jobscript arguments
    parser.add_argument("-jsched","--jsched", help="job scheduling system. Choices: SLURM or SGE",action="store_true")
    parser.add_argument("-jname","--jname", help="jobs main identifier",action="store_true")
    parser.add_argument("-memory","--memory", help="memory reserved per thread for job file in G(default: 2G)e.g.2",action="store_true")
    parser.add_argument("-wtime","--wtime", help="wall time requested in hours for queueing system (default: 168hrs) e.g. 8",action="store_true")
    parser.add_argument("-queue","--queue", help="queue name e.g gpus",action="store_true")
    parser.add_argument("-gpus","--gpus", help="number of GPUS (default: 1)",action="store_true")
    parser.add_argument("-cpus","--cpus", help="number of CPUs (default: 1)",action="store_true")
    parser.add_argument("-modules","--modules", help="modules to be loaded for the calculation",action="store_true")
    parser.add_argument("-joption","--joption", help="additional options for jobscript",action="store_true")
    parser.add_argument("-jcommand","--jcommand", help="additional commands for jobscript",action="store_true")
    # database search arguments
    parser.add_argument("-dbsim","--dbsim", help="SMILES/ligand/file for similarity search",action="store_true")
    parser.add_argument("-dbcatoms","--dbcatoms", help="connection atoms for similarity search",action="store_true")
    parser.add_argument("-dbresults","--dbresults", help="how many results for similary search or screening",action="store_true")
    parser.add_argument("-dboutputf","--dboutputf", help="output file for search results",action="store_true")
    parser.add_argument("-dbbase","--dbbase", help="database for search",action="store_true")
    parser.add_argument("-dbsmarts","--dbsmarts", help="SMARTS string for screening",action="store_true")
    parser.add_argument("-dbfinger","--dbfinger", help="fingerprint for similarity search",action="store_true")
    parser.add_argument("-dbatoms","--dbatoms", help="number of atoms to be used in screening",action="store_true")
    parser.add_argument("-dbbonds","--dbbonds", help="number of bonds to be used in screening",action="store_true")
    parser.add_argument("-dbarbonds","--dbarbonds", help="Number of aromatic bonds to be used in screening",action="store_true")
    parser.add_argument("-dbsbonds","--dbsbonds", help="number of single bonds to be used in screening",action="store_true")
    parser.add_argument("-dbmw","--dbmw", help="molecular weight to be used in screening",action="store_true")
    # post-processing arguments
    parser.add_argument("-postp","--postp",help="post process results",action="store_true")
    parser.add_argument("-postqc","--postqc", help="quantum chemistry code used. Choices: TeraChem or GAMESS",action="store_true") 
    parser.add_argument("-postdir","--postdir", help="directory with results",action="store_true") 
    parser.add_argument("-pres","--pres",help="generate calculations summary",action="store_true")
    parser.add_argument("-pdeninfo","--pdeninfo",help="calculate average properties for electron density",action="store_true")
    parser.add_argument("-pcharge","--pcharge",help="calculate charges",action="store_true")
    parser.add_argument("-pgencubes","--pgencubes",help="generate cubefiles",action="store_true")
    parser.add_argument("-pwfninfo","--pwfninfo",help="get information about wavefunction",action="store_true")
    parser.add_argument("-pdeloc","--pdeloc",help="get delocalization and localization indices",action="store_true")
    parser.add_argument("-porbinfo","--porbinfo",help="get information about MO",action="store_true")
    parser.add_argument("-pnbo","--pnbo",help="post process nbo analysis",action="store_true")
    #parser.add_argument("-pdorbs","--pdorbs",help="get information on metal d orbitals",action="store_true")
    # auxiliary
    parser.add_argument("-dbsearch","--dbsearch",action="store_true") # flag for db search
    parser.add_argument("-checkdirt","--checkdirt",action="store_true") # directory removal check flag
    parser.add_argument("-checkdirb","--checkdirb",action="store_true") # directory removal check flag 2
    parser.add_argument("-jid","--jid",action="store_true")           # job id for folders
    parser.add_argument("-gui","--gui",action="store_true")           # gui placeholder
    # calculations summary (terapostp gampost)
    # nbo
    # charges
    # wavefunction - cube files
    # deloc indices - basin analysis
    # moldparse
    args=parser.parse_args()
    return args
