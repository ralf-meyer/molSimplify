# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

#########################################################
############ This script generates input  ###############
##########  for Quantum Chem calculations   #############
#########################################################
from Classes.globalvars import *
import glob, sys, os, shutil
from Classes.mol3D import mol3D
from Classes.atom3D import atom3D
from Classes.globalvars import *

#####################################################
########## This module generates input  #############
##########  for TeraChem calculations   #############
#####################################################
### generate multiple runs if multiple methods requested ###
def multitcgen(args,strfiles):
    jobdirs = []
    method = False
    if args.method and len(args.method) > 1:
        methods = args.method
        for method in methods:
            jobdirs.append(tcgen(args,strfiles,method))
    else:
        jobdirs.append(tcgen(args,strfiles,method))
    # remove original files
    for xyzf in strfiles:
        os.remove(xyzf+'.xyz')
        os.remove(xyzf+'.molinp')
    return jobdirs

### generate terachem input files ###
def tcgen(args,strfiles,method):
    # global variables
    globs = globalvars()
    jobdirs = []
    coordfs = []
    # Initialize the jobparams dictionary with mandatory/useful keywords.
    jobparams={'run': 'minimize',
           'timings': 'yes',
           'min_coordinates':'cartesian',
           'maxit': '500',
           'scrdir': './scr',
           'method': 'b3lyp',
           'basis': 'lacvps_ecp',
           'spinmult': '1',
           'charge': '0',
           'gpus': '1',
            }
    # if multiple methods requested generate c directories
    # Overwrite plus add any new dictionary keys from commandline input.       
    for xyzf in strfiles:
        rdir = xyzf.rsplit('/',1)[0]
        xyzft = xyzf.rsplit('/',1)[-1]
        xyzf += '.xyz'
        coordfs.append(xyzf.rsplit('/',1)[-1])
        coordname = xyzft
        # Setting jobname for files + truncated name for queue.
        if len(coordname) > 10:
            nametrunc=coordname[0:6]+coordname[-4:]
        else:
            nametrunc=coordname
        if not os.path.exists(rdir+'/'+nametrunc):
            os.mkdir(rdir+'/'+nametrunc) 
        mdir = rdir+'/'+nametrunc
        if method:
            if method[0]=='U' or method[0]=='u':
                mmd = '/'+method[1:]
            else:
                mmd = '/'+method
            mdir = rdir+'/'+nametrunc+mmd
            if not os.path.exists(mdir):
                os.mkdir(mdir)
        jobdirs.append(mdir)
        shutil.copy2(xyzf,mdir)
        shutil.copy2(xyzf.replace('.xyz','.molinp'),mdir.replace('.xyz','.molinp'))
    # parse extra arguments
    # Method parsing, does not check if a garbage method is used here:
    unrestricted=False
    if method:
        jobparams['method'] = method
        if ('u' or 'U') in method[0]:
            # Unrestricted calculation
            unrestricted=True
        else:
            # Restricted calculation
            unrestricted=False
            if args.spin and int(args.spin) > 1:
                jobparams['method'] = 'u'+method
                unrestricted=True
    else:
        if args.spin and int(args.spin) > 1:
                jobparams['method'] = 'ub3lyp'
                unrestricted=True
        else:
            jobparams['method'] = 'b3lyp'
    if (args.runtyp and 'energy' in args.runtyp.lower()):
        jobparams['run'] = 'energy'
    elif (args.runtyp and 'ts' in args.runtyp.lower()):
        jobparams['run'] = 'ts'
    elif (args.runtyp and 'gradient' in args.runtyp.lower()):
        jobparams['run'] = 'gradient'
    if (args.gpus):
        jobparams['gpus'] = args.gpus
    if (args.dispersion):
        jobparams['dispersion']=args.dispersion
    # Just carry over spin and charge keywords if they're set. Could do checks, none for now.
    if args.spin:
        jobparams['spinmult']=args.spin
    if args.charge:
        if args.bcharge:
            args.charge = int(args.charge)+int(args.bcharge)
        jobparams['charge']=args.charge
    # Check for existence of basis and sanitize name
    if args.basis:
        ecp=False # Flag not currently used, for deciding gpus_ecp code or not later. Can always specify with 'extra' command
        if '*' in args.basis:
            jobparams['basis']=args.basis.replace('*','s')
        else:
            jobparams['basis']=args.basis
    # Overwrite plus add any new dictionary keys from commandline input.       
    if args.qoption:
        if len(args.qoption)%2!=0:
                print 'WARNING: wrong number of arguments in -qoption'
        else:
            for elem in range(0,int(0.5*len(args.qoption))):
                key,val=args.qoption[2*elem],args.qoption[2*elem+1]
                jobparams[key]=val
    # Extra keywords for unrestricted. 
    if unrestricted:
       # If running unrestricted, assume convergence will be more difficult for now.
       jobparams['scf']='diis+a' 
       if not jobparams.has_key('levelshift'):
          jobparams['levelshift']='yes'
       elif jobparams['levelshift'] != 'yes':
          print("Warning! You're doing an unrestricted calculation but have set levelshift = %s" %(jobparams['levelshift']))
       if not jobparams.has_key('levelshiftvala'):
          jobparams['levelshiftvala']='1.0'
       if not jobparams.has_key('levelshiftvalb'):
          jobparams['levelshiftvalb']='0.1'
    # Now we're ready to start building the input file
    for i,jobd in enumerate(jobdirs):
        output=open(jobd+'/terachem_input','w')
        output.write('# file created with %s\n' % globs.PROGRAM)
        jobparams['coordinates'] = coordfs[i]
        for keys in jobparams.keys():
            output.write('%s %s\n' %(keys,jobparams[keys]))
        output.write('end\n')
        output.close()
    return jobdirs

#####################################################
########## This module generates input  #############
##########    for GAMESS calculations   #############
#####################################################

### convert xyz to gxyz ###
def xyz2gxyz(filename):
    # convert normal xyz file to gamess format
    mol=mol3D() # create mol3D object
    mol.readfromxyz(filename) # read molecule
    gfilename = filename.replace('.xyz','.gxyz') # new file name
    mol.writegxyz(gfilename) # write gamess formatted xyz file
    return gfilename.split('.gxyz')[0]

### generate multiple runs if multiple methods requested ###
def multigamgen(args,strfiles):
    method = False
    jobdirs=[]
    if args.method and len(args.method) > 1:
        methods = args.method
        for method in methods:
            jobdirs.append(gamgen(args,strfiles,method))
    else:
        jobdirs.append(gamgen(args,strfiles,method))
    # remove original files
    for xyzf in strfiles:
        os.remove(xyzf+'.xyz')
        os.remove(xyzf+'.gxyz')
        os.remove(xyzf+'.molinp')
    return jobdirs

### generate input files for gamess###
def gamgen(args,strfiles,method):
    # get global variables
    globs = globalvars()
    jobdirs = []
    coordfs = []
    # Initialize the jobparams dictionary with mandatory/useful keywords.
    jobparams={'RUNTYP': 'OPTIMIZE',
           'GBASIS': 'N21',
           'MAXIT': '500',
           'DFTTYP': 'B3LYP',
           'SCFTYP': 'UHF',
           'ICHARG': '0',
           'MULT': '1',
            }
    # Overwrite plus add any new dictionary keys from commandline input.       
    for xyzf in strfiles:
        # convert to "gamess format"
        xyzf = xyz2gxyz(xyzf+'.xyz')
        rdir = xyzf.rsplit('/',1)[0]
        xyzft = xyzf.rsplit('/',1)[-1]
        xyzf += '.gxyz'
        coordfs.append(xyzf)
        coordname = xyzft
        # Setting jobname for files + truncated name for queue.
        if len(coordname) > 10:
            nametrunc=coordname[0:6]+coordname[-4:]
        else:
            nametrunc=coordname
        if not os.path.exists(rdir+'/'+nametrunc):
            os.mkdir(rdir+'/'+nametrunc) 
        mdir = rdir+'/'+nametrunc
        if method:
            if method[0]=='U' or method[0]=='u':
                mmd = '/'+method[1:]
            else:
                mmd = '/'+method
                jobparams['SCFTYP']='RHF'
            mdir = rdir+'/'+nametrunc+mmd
            if not os.path.exists(mdir):
                os.mkdir(mdir)
        jobdirs.append(mdir)
        shutil.copy2(xyzf,mdir)
        shutil.copy2(xyzf.replace('.gxyz','.molinp'),mdir.replace('.gxyz','.molinp'))
    if method:
        if method[0]=='U' or method[0]=='u':
            method = method[1:]
    # Just carry over spin and charge keywords if they're set. Could do checks, none for now.
    if args.spin:
       jobparams['MULT']=str(args.spin)
    if args.charge:
       jobparams['ICHARG']=str(args.charge)
    # Check for existence of basis and sanitize name
    if args.gbasis:
          jobparams['GBASIS']=args.gbasis.upper()
    if args.ngauss:
          jobparams['NGAUSS']=args.ngauss.upper()
    if method:
          jobparams['DFTTYP']=method.upper()
    if (args.runtyp and 'en' in args.runtyp.lower()):
        jobparams['run'] = 'ENERGY'
    elif (args.runtyp and 'ts' in args.runtyp.lower()):
        jobparams['run'] = 'SADPOINT'
    # Now we're ready to start building the input file and the job script
    for i,jobd in enumerate(jobdirs):
        output=open(jobd+'/gam.inp','w')
        f=open(coordfs[i])
        s = f.read() # read coordinates
        f.close()
        jobparams['coordinates'] = s
        output.write('! File created using %s\n' % globs.PROGRAM)
        # write $BASIS block
        output.write(' $BASIS ')
        if args.ngauss:
            output.write(' GBASIS='+jobparams['GBASIS'])
            output.write(' NGAUSS='+jobparams['NGAUSS'])
        else:
            output.write(' GBASIS='+jobparams['GBASIS'])
        if args.ndfunc:
            output.write(' NDFUNC='+args.ndfunc)
        if args.npfunc:
            output.write(' NPFUNC='+args.npfunc)
        output.write(' $END\n')
        # write $SYSTEM block
        output.write(' $SYSTEM ')
        # check if MWORDS specified by the user
        if not args.sysoption or not ('MWORDS' in args.sysoption):
            output.write(' MWORDS=16')
        # write additional options
        if (args.sysoption):
            if len(args.sysoption)%2 > 0:
                print 'WARNING: wrong number of arguments in -sysoption'
            else:
                for elem in range(0,int(0.5*len(args.sysoption))):
                    key,val=args.sysoption[2*elem],args.sysoption[2*elem+1]
                    output.write(' '+key+'='+val+' ')
        output.write(' $END\n')
        # write CONTRL block
        output.write(' $CONTRL SCFTYP='+jobparams['SCFTYP']+' DFTTYP=')
        output.write(jobparams['DFTTYP']+' RUNTYP='+jobparams['RUNTYP'])
        output.write('\n  ICHARG='+jobparams['ICHARG']+' MULT=')
        # check if CC basis set specified and add spherical
        if 'CC' in jobparams['GBASIS']:
            output.write(jobparams['MULT']+' ISPHER=1\n')
        else:
            output.write(jobparams['MULT']+'\n')
        # write additional options
        if (args.ctrloption):
            if len(args.ctrloption)%2 > 0:
                print 'WARNING: wrong number of arguments in -ctrloption'
            else:
                for elem in range(0,int(0.5*len(args.ctrloption))):
                    key,val=args.ctrloption[2*elem],args.ctrloption[2*elem+1]
                    output.write(' '+key+'='+val+' ')
        output.write(' $END\n')
        # write $SCF block
        output.write(' $SCF ')
        # check if options specified by the user
        if not args.scfoption or not ('DIRSCF' in args.scfoption):
            output.write(' DIRSCF=.TRUE.')
        if not args.scfoption or not ('DIIS' in args.scfoption):
            output.write(' DIIS=.TRUE.')
        if not args.scfoption or not ('SHIFT' in args.scfoption):
                output.write(' SHIFT=.TRUE.')
        # write additional options
        if (args.scfoption):
            if len(args.scfoption)%2!=0:
                print 'WARNING: wrong number of arguments in -scfoption'
            else:
                for elem in range(0,int(0.5*len(args.scfoption))):
                    key,val=args.scfoption[2*elem],args.scfoption[2*elem+1]
                    output.write(' '+key+'='+val+' ')
        output.write(' $END\n')
        # write $STATPT block
        output.write(' $STATPT ')
        # check if NSTEP specified by the user
        if not args.statoption or not ('NSTEP' in args.statoption):
            output.write(' NSTEP=100')
        # write additional options
        if (args.statoption):
            if len(args.statoption)%2 > 0:
                print 'WARNING: wrong number of arguments in -statoption'
            else:
                for elem in range(0,int(0.5*len(args.statoption))):
                    key,val=args.statoption[2*elem],args.statoption[2*elem+1]
                    output.write(' '+key+'='+val+' ')
        output.write(' $END\n')
        # write $DATA block
        output.write(' $DATA\n')
        output.write(jobparams['coordinates']+' $END\n')
        output.close()
    return jobdirs

    
####################################################
########## This module generates input  ############
##########    for QChem calculations   #############
####################################################
### generate multiple runs if multiple methods requested ###
def multiqgen(args,strfiles):
    method = False
    jobdirs=[]
    if args.method and len(args.method) > 1:
        methods = args.exchange
        for method in methods:
            jobdirs.append(qen(args,strfiles,method))
    else:
        jobdirs.append(qgen(args,strfiles,method))
    # remove original files
    for xyzf in strfiles:
        os.remove(xyzf+'.xyz')
        os.remove(xyzf+'.molinp')
    return jobdirs

### generate input files for qchem ###
def qgen(args,strfiles,method):
    # get global variables
    globs = globalvars()
    jobdirs = []
    coordfs = []
    # Initialize the jobparams dictionary with mandatory/useful keywords.
    jobparams={'UNRESTRICTED': 'true',
           'BASIS': 'lanl2dz',
           'JOBTYPE': 'opt',
            'EXCHANGE': 'b3lyp',
           'CORRELATION': 'none',
           'MAX_SCF_CYCLES': '500',
           'GEOM_OPT_MAX_CYCLES': '1000',
           'SYMMETRY': 'off',
           'PRINT_ORBITALS': 'true',
           'CHARGE':'1',
           'SPIN':'1',
            }
    # Overwrite plus add any new dictionary keys from commandline input.       
    for xyzf in strfiles:
        rdir = xyzf.rsplit('/',1)[0]
        xyzft = xyzf.rsplit('/',1)[-1]
        xyzf += '.xyz'
        coordfs.append(xyzf.rsplit('/',1)[-1])
        coordname = xyzft
        # Setting jobname for files + truncated name for queue.
        if len(coordname) > 10:
            nametrunc=coordname[0:6]+coordname[-4:]
        else:
            nametrunc=coordname
        if not os.path.exists(rdir+'/'+nametrunc):
            os.mkdir(rdir+'/'+nametrunc) 
        mdir = rdir+'/'+nametrunc
        if method:
            mmd = '/'+method
            mdir = rdir+'/'+nametrunc+mmd
            if not os.path.exists(mdir):
                os.mkdir(mdir)
        jobdirs.append(mdir)
        shutil.copy2(xyzf,mdir)
        shutil.copy2(xyzf.replace('.xyz','.molinp'),mdir.replace('.xyz','.molinp'))
    # Check for existence of basis and sanitize name
    if args.basis and len(args.basis) > 1:
        jobparams['BASIS']=args.basis
    if args.correlation and len(args.correlation) > 1:
        jobparams['CORRELATION']=args.correlation
    if method and len(method) > 1:
        jobparams['EXCHANGE']=method
    if not args.unrestricted:
        jobparams['UNRESTRICTED']='false'
    if (args.runtyp and 'en' in args.runtyp.lower()):
        jobparams['run'] = 'SP'
    elif (args.runtyp and 'ts' in args.runtyp.lower()):
        jobparams['run'] = 'TS'
    # Just carry over spin and charge keywords if they're set. Could do checks, none for now.
    if args.spin:
       jobparams['SPIN']=args.spin
    if args.charge:
       jobparams['CHARGE']=args.charge
    # Now we're ready to start building the input file and the job script
    for i,jobd in enumerate(jobdirs):
        output=open(jobd+'/qch.inp','w')
        f=open(jobd+'/'+coordfs[i])
        s0 = f.readlines()[2:] # read coordinates
        f.close()
        # if separate split to two molecules
        if args.bsep and '--' in ''.join(s0):
            idxsplit = [isdx for isdx, ss in enumerate(s0) if '--' in ss][0]
            s = '--\n'+jobparams['CHARGE']+' '+jobparams['SPIN']+'\n'
            s += ''.join(s0[:idxsplit])
            s += '--\n0 1\n'
            s += ''.join(s0[idxsplit+3:])
        else:
            s = s0
        # write rem block
        output.write('$rem\nUNRESTRICTED\t\t' + jobparams['UNRESTRICTED'])
        output.write('\nBASIS\t\t'+jobparams['BASIS']+'\nJOBTYPE\t\t'+jobparams['JOBTYPE'])
        output.write('\nEXCHANGE\t\t'+jobparams['EXCHANGE']+'\nCORRELATION\t\t')
        output.write(jobparams['CORRELATION']+'\nMAX_SCF_CYCLES\t\t')
        output.write(jobparams['MAX_SCF_CYCLES']+'\nGEOM_OPT_MAX_CYCLES\t\t')
        output.write(jobparams['GEOM_OPT_MAX_CYCLES']+'\nSYMMETRY\t\t'+jobparams['SYMMETRY'])
        output.write('\nPRINT_ORBITALS\t\t'+jobparams['PRINT_ORBITALS']+'\n')
        # write additional options
        if (args.remoption):
            if len(args.remoption)%2 > 0:
                print 'WARNING: wrong number of arguments in -remoption'
            else:
                for elem in range(0,int(0.5*len(args.remoption))):
                    key,val=args.remoption[2*elem],args.remoption[2*elem+1]
                    output.write(key+'\t\t'+val+'\n')
        output.write('$end\n\n')
        # write $molecule block
        output.write('$molecule\n'+jobparams['CHARGE']+' '+jobparams['SPIN']+'\n')
        output.write(''.join(s)+'$end')
        output.close()
    return jobdirs
