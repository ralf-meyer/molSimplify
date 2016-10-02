# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

#################################################################
######## This scripts adds new molecules to our database ########
#################################################################

# import custom modules
from geometry import *
from io import *
from Classes.globalvars import *
# import std modules
import os, sys, subprocess, re, unicodedata
import pybel, openbabel, random, shutil


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
    licores = readdict(globs.installdir+'/Ligands/ligands.dict')
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
        lig,emsg = lig_load(globs.installdir+'/',smimol,licores)
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
            shutil.copy2(smimol,globs.installdir+'/Ligands/'+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+shortname+','+css+','+grp+','+ffopt
        elif '.xyz' in smimol:
            shutil.copy2(smimol,globs.installdir+'/Ligands/'+sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css+','+grp+','+ffopt
        elif lig.OBmol:
            # write smiles file in Ligands directory
            lig.OBmol.write('smi',globs.installdir+'/Ligands/'+sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','+css+','+grp+','+ffopt
        else:
            # write xyz file in Ligands directory
            lig.writexyz(globs.installdir+'/Ligands/'+sminame+'.xyz') # write xyz file
        # update dictionary
        f = open(globs.installdir+'/Ligands/ligands.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Ligands/ligands.dict','w')
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
    mcores = readdict(globs.installdir+'/Cores/cores.dict')
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
        core,emsg = core_load(globs.installdir+'/',smimol,mcores)
        if emsg:
            return emsg
        core.convert2mol3D() # convert to mol3D
        # write xyz file in Cores directory
        # new entry for dictionary
        if '.mol' in smimol:
            shutil.copy2(smimol,globs.installdir+'/Cores/'+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+css+','+'1'
        elif '.xyz' in smimol:
            shutil.copy2(smimol,globs.installdir+'/Cores/'+sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        else:
            core.writexyz(globs.installdir+'/Cores/'+sminame+'.xyz') # write xyz file
            # new entry for dictionary
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        # update dictionary
        f = open(globs.installdir+'/Cores/cores.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Cores/cores.dict','w')
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
    bindcores = readdict(globs.installdir+'/Bind/bind.dict')
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
        bind,bsmi,emsg = bind_load(globs.installdir+'/',smimol,bindcores)
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
            shutil.copy2(smimol,globs.installdir+'/Bind/'+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+shortname+','+css
        elif '.xyz' in smimol:
            shutil.copy2(smimol,globs.installdir+'/Bind/'+sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css
        elif bind.OBmol:
            # write smiles file in Bind species directory
            bind.OBmol.write('smi',globs.installdir+'/Bind/'+sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','+css
        else:
            # write xyz file in Bind species directory
            bind.writexyz(globs.installdir+'/Bind/'+sminame+'.xyz') # write xyz file
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css
        # update dictionary
        f = open(globs.installdir+'/Bind/bind.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Bind/bind.dict','w')
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
    # convert to unicode
    sminame = unicodedata.normalize('NFKD',sminame).encode('ascii','ignore')
    if ropt==1:
        # update dictionary
        f = open(globs.installdir+'/Ligands/ligands.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Ligands/ligands.dict','w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(globs.installdir+'/Ligands/'+sss[1].split(',')[0])
        f.close()
    elif ropt==0:
        mcores = readdict(globs.installdir+'/Cores/cores.dict')
        # update dictionary
        f = open(globs.installdir+'/Cores/cores.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Cores/cores.dict','w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(globs.installdir+'/Cores/'+sss[1].split(',')[0])
        f.close()
    elif ropt==2:
        bindcores = readdict(globs.installdir+'/Bind/bind.dict')
        # update dictionary
        f = open(globs.installdir+'/Bind/bind.dict','r')
        ss = f.read().splitlines()
        f.close()
        f = open(globs.installdir+'/Bind/bind.dict','w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(globs.installdir+'/Bind/'+sss[1].split(',')[0])
        f.close()
    return emsg
