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
import os, sys, subprocess, re, unicodedata
import pybel, openbabel, random, shutil
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
    licores = readdict(resource_filename(Requirement.parse("molSimplify"),"molSimplify/Ligands/ligands.dict"))
    ligands_folder = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Ligands")
    print("ligands_folder is : " + str(ligands_folder))
   # licores = readdict(globs.installdir+'/Ligands/ligands.dict')
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
            shutil.copy2(smimol,ligans_folder + sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css+','+grp+','+ffopt
        elif lig.OBmol: 
            # write smiles file in Ligands directory
            lig.OBmol.write('smi',ligands_folder + sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','+css+','+grp+','+ffopt
        else:
            # write xyz file in Ligands directory
            lig.writexyz(ligands_folder+sminame+'.xyz') # write xyz file
        # update dictionary
        lipath = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Ligands/ligands.dict")


#        f = open(globs.installdir+'/Ligands/ligands.dict','r')
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
    mcores = readdict(resource_filename(Requirement.parse("molSimplify"),"molSimplify/Cores/cores.dict"))
    core_dir  = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Cores")
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
            shutil.copy2(smimol,core_dir+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+css+','+'1'
        elif '.xyz' in smimol:
            shutil.copy2(smimol,core_dir + sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        else:
            core.writexyz(core_dir +sminame+'.xyz') # write xyz file
            # new entry for dictionary
            snew = sminame+':'+sminame+'.xyz,'+css+','+'1'
        # update dictionary
        cpath = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Cores/cores.dict")
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
    bindcores = readdicct(resource_filename(Requirement.parse("molSimplify"),"molSimplify/Bind/bind.dict"))
    bind_dir = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Bind/bind.dict")
    print(bind_dir)

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
            shutil.copy2(smimol,bind_dir+sminame+'.mol')
            snew = sminame+':'+sminame+'.mol,'+shortname+','+css
        elif '.xyz' in smimol:
            shutil.copy2(smimol,bind_dir +sminame+'.xyz')
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css
        elif bind.OBmol:
            # write smiles file in Bind species directory
            bind.OBmol.write('smi',bind_dir +sminame+'.smi')
            snew = sminame+':'+sminame+'.smi,'+shortname+','+css
        else:
            # write xyz file in Bind species directory
            bind.writexyz(bind_dir +sminame+'.xyz') # write xyz file
            snew = sminame+':'+sminame+'.xyz,'+shortname+','+css
        # update dictionary
        bpath = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Bind/bind.dict")
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
    # convert to unicode
    sminame = unicodedata.normalize('NFKD',sminame).encode('ascii','ignore')
    lipath = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Ligands/ligands.dict")
    li_dir = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Ligands")
    core_path = (resource_filename(Requirement.parse("molSimplify"),"molSimplify/Cores/cores.dict"))
    core_dir  = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Cores")
    bind_path = (resource_filename(Requirement.parse("molSimplify"),"molSimplify/Bind/bind.dict"))
    bind_dir  = resource_filename(Requirement.parse("molSimplify"),"molSimplify/Bind")


    if ropt==1:
        # update dictionary
        f = open(lipath,'r')
        ss = f.read().splitlines()
        f.close()
        f = open(lipath,'w')
        ssort = sorted(ss[1:])
        f.write(ss[0]+'\n')
        for s in ssort:
            sss = s.split(':')
            if sminame!=sss[0]:
                f.write(s+'\n')
            else:
                os.remove(li_dir + sss[1].split(',')[0])
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
                os.remove(core_dir+sss[1].split(',')[0])
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
                os.remove(bind_dir+sss[1].split(',')[0])
        f.close()
    return emsg
