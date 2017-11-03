# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

#################################################################
######## This scripts adds new molecules to our database ########
#################################################################

# import custom modules
from geometry import *
from io import *
from molSimplify.Classes.globalvars import *
# import std modules
from molSimplify.Classes.mWidgets import *

import os, sys, subprocess, re, unicodedata
import openbabel, random, shutil
from pkg_resources import resource_filename, Requirement

###############################
### adds to ligand database ###
###############################
def addtoldb(smimol,sminame,smident,smicat,smigrps,smictg,ffopt):
    #  INPUT
    #   - smimol: SMILES string or molecule file to be added
    #   - sminame: name of ligand for key in dictionary
    #   - smident: denticity of the ligand
    #   - smicat: connection atoms
    #  OUTPUT
    #   - emsg: error messages
    emsg = False
    globs = globalvars()
    if not globs.custom_path  or not os.path.exists(str(globs.custom_path)):
	## here, we need to give a path, try and
	## get one from cmd line:
#        if args.gui:
	    ## the GUI is handled in the GUI script,
	    ## this should not reached during normal
#	    ## operation (it's n
#            qqb = mQDialogWarn('Warning','No custom user-writable directory found, addition not  not possible.')
#            qqb.setParent(args.gui.DBWindow)
#            raise Exception('Custom path not set!')
#	else: 
	print('To add to database, you need to set a custom path. Please enter a writeable file path:')
	new_path = input('path=')
	globs.add_custom_path(new_path)
       	copy_to_custom_path()

    lipath = globs.custom_path + "/Ligands/ligands.dict"
    licores = readdict(lipath)
    ligands_folder = globs.custom_path + "/Ligands/"
    print("ligands_folder is : " + str(ligands_folder))
    # check if ligand exists
    if sminame in licores.keys():
        emsg = 'Ligand '+sminame+' already existing in ligands database.' 
        emsg += ' To replace, delete the existing entry first.'
        return emsg
    else:
        # get connection atoms
        ccats = filter(None,re.split(' |,|\t',smicat))
        # get groups
        groups = filter(None,re.split(' |,|\t',smigrps))
        grp = 'all '+' '.join(groups)
        grp += ' '+smictg
        if smicat=='':
            cats = range(0,int(smident))
        else:
            cats = [int(a)-1 for a in ccats]
        cs = [str(a) for a in cats]
        css = ' '.join(cs)
        # convert to unicode
        smimol = unicodedata.normalize('NFKD',smimol).encode('ascii','ignore')
        sminame = unicodedata.normalize('NFKD',sminame).encode('ascii','ignore')
        if '~' in smimol:
            smimol = smimol.replace('~',os.expanduser('~'))
        # convert ligand from smiles/file
        lig,emsg = lig_load(smimol,licores)
        if emsg:
            return emsg
        lig.convert2mol3D() # convert to mol3D
        # create shortname
        if len(sminame) > 5:
            shortname = sminame[0:3]+sminame[-2:]
        else:
            shortname = sminame
        # new entry for dictionary
        if '.mol' in smimol:
            shutil.copy2(smimol,ligands_folder + sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+shortname+','+css+','+grp+','+ffopt
        elif '.xyz' in smimol:
            shutil.copy2(smimol,ligands_folder + sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css+','+grp+','+ffopt
        elif lig.OBmol: 
            # write smiles file in Ligands directory
            lig.OBmol.write('smi',ligands_folder + sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','+css+','+grp+','+ffopt
        else:
            # write xyz file in Ligands directory
            lig.writexyz(ligands_folder+sminame+'.xyz') # write xyz file
        # update dictionary


        f = open(lipath,'r')

        ss = f.read().splitlines()
        f.close()
        f = open(lipath,'w')
        ss.append(snew)
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            f.write(s+'\n')
        f.close()
    return emsg
    
##############################
### adds to cores database ###
##############################
def addtocdb(smimol,sminame,smicat):
    #  INPUT
    #   - smimol: SMILES string or molecule file to be added
    #   - sminame: name of core for key in dictionary
    #   - smicat: connection atoms
    emsg = False
    globs = globalvars()
    if not globs.custom_path  or not os.path.exists(str(globs.custom_path)):
    	## here, we need to give a path, try and
    	## get one from cmd line:
    	## the GUI is handled in the GUI script,
    	## this should not reached during normal
    	## operation 
    	print('To add to database, you need to set a custom path. Please enter a writeable file path:')
    	new_path = input('path=')
    	globs.add_custom_path(new_path)
    	copy_to_custom_path()

    cpath = globs.custom_path + "/Cores/cores.dict"
    mcores = readdict(cpath)
    cores_folder = globs.custom_path + "/Cores/"
    # check if core exists
    if sminame in mcores.keys():
        emsg = 'Core '+sminame+' already existing in core database.'
        return emsg
    else:
        # get connection atoms
        ccats = filter(None,re.split(' |,|\t',smicat))
        cats = [int(a)-1 for a in ccats]
        if len(cats)==0:
            cats=[0]
        cs = [str(a) for a in cats]
        css = ' '.join(cs)
        # convert to unicode
        smimol = unicodedata.normalize('NFKD',smimol).encode('ascii','ignore')
        if '~' in smimol:
            smimol = smimol.replace('~',os.expanduser('~'))
        # convert ligand from smiles/file
        core,emsg = core_load(smimol,mcores)
        if emsg:
            return emsg
        core.convert2mol3D() # convert to mol3D
        # write xyz file in Cores directory
        # new entry for dictionary
        if '.mol' in smimol:
            shutil.copy2(smimol,cores_folder+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+css+','+'1'
        elif '.xyz' in smimol:
            shutil.copy2(smimol,cores_folder + sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        else:
            core.writexyz(cores_folder +sminame+'.xyz') # write xyz file
            # new entry for dictionary
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        # update dictionary
        f = open(cpath,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(cpath,'w')
        ss.append(snew)
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            f.write(s+'\n')
        f.close()
    return emsg
    
########################################
### adds to binding species database ###
########################################
def addtobdb(smimol,sminame):
    #  INPUT
    #   - smimol: SMILES string or molecule file to be added
    #   - sminame: name of binding species for key in dictionary
    globs = globalvars()
    if not globs.custom_path  or not os.path.exists(str(globs.custom_path)):
    	## here, we need to give a path, try and
    	## get one from cmd line:
    	## the GUI is handled in the GUI script,
    	## this should not reached during normal
    	## operation 
    	print('To add to database, you need to set a custom path. Please enter a writeable file path:')
    	new_path = input('path=')
    	globs.add_custom_path(new_path)
    	copy_to_custom_path()

    bpath = globs.custom_path + "/Bind/bind.dict"
    bindcores = readdict(bpath)
    bind_folder = globs.custom_path + "/Bind/"


    # check if binding species exists
    if sminame in bindcores.keys():
        emsg = 'Molecule '+sminame+' already existing in binding species database.'
        return emsg
    else:
        # convert to unicode
        smimol = unicodedata.normalize('NFKD',smimol).encode('ascii','ignore')
        sminame = unicodedata.normalize('NFKD',sminame).encode('ascii','ignore')
        if '~' in smimol:
            smimol = smimol.replace('~',os.expanduser('~'))
        # convert ligand from smiles/file
        bind,bsmi,emsg = bind_load(smimol,bindcores)
        if emsg:
            return emsg
        bind.convert2mol3D() # convert to mol3D
                # new entry for dictionary
                # create shortname
        if len(sminame) > 5:
            shortname = sminame[0:3]+sminame[-2:]
        else:
            shortname = sminame
        if '.mol' in smimol:
            shutil.copy2(smimol,bind_folder+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+shortname+','
        elif '.xyz' in smimol:
            shutil.copy2(smimol,bind_folder +sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','
        elif bind.OBmol:
            # write smiles file in Bind species directory
            bind.OBmol.write('smi',bind_folder +sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','
        else:
            # write xyz file in Bind species directory
            bind.writexyz(bind_folder +sminame+'.xyz') # write xyz file
            snew = sminame+':'+sminame+'.xyz,'+shortname+','
        # update dictionary
        f = open(bpath,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(bpath,'w')
        ss.append(snew)
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            f.write(s+'\n')
        f.close()
    return emsg
    
############################
### remove from database ###
############################
def removefromDB(sminame,ropt):
    #  INPUT
    #   - sminame: name of molecule for key in dictionary
    #  OUTPUT
    #   - emsg: error messages
    emsg = False
    globs = globalvars()
    if not globs.custom_path  or not os.path.exists(str(globs.custom_path)):
    	## here, we need to give a path, try and
    	## get one from cmd line:
    	## the GUI is handled in the GUI script,
    	## this should not reached during normal
    	## operation 
    	print('To database, you need to set a custom path. Please enter a writeable file path:')
    	new_path = input('path=')
    	globs.add_custom_path(new_path)
    	copy_to_custom_path()
    li_path =  globs.custom_path + "/Ligands/ligands.dict"
    li_folder =  globs.custom_path +  "/Ligands/"
    core_path =  globs.custom_path +  "/Cores/cores.dict"
    core_dir  =  globs.custom_path +  "/Cores/"
    bind_path =  globs.custom_path + "/Bind/bind.dict"
    bind_folder  = globs.custom_path +"/Bind/"

    # convert to unicode
    sminame = unicodedata.normalize('NFKD',sminame).encode('ascii','ignore')

    if ropt==1:
        # update dictionary
        f = open(li_path,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(li_path,'w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(li_folder + sss[1].split(',')[0])
        f.close()
    elif ropt==0:
        mcores = readdict(core_path)
        # update dictionary
        f = open(core_path,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(core_path,'w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(core_folder+sss[1].split(',')[0])
        f.close()
    elif ropt==2:
        bindcores = readdict(bind_path)
        # update dictionary
        f = open(bind_path,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(bind_path,'w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(bind_folder+sss[1].split(',')[0])
        f.close()
    return emsg
