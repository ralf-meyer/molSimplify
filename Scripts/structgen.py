# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## This script generates a collection  #############
#####  of randomly placed binding species around   #######
########   a functionalized ferrocene core   #############
######## and then creates the required input #############
######## and job files for running terachem  #############
########  calculations on these structures   #############
##########################################################

# import custom modules
from geometry import *
from io import *
from Classes.globalvars import *
# import standard modules
import os, sys
import pybel, openbabel, random, itertools
from numpy import log, arccos, cross, dot, pi

######################
### Euclidean norm ###
######################
def norm(u):
    # INPUT
    #   - u: n-element list
    # OUTPUT
    #   - d: Euclidean norm
    d = 0.0
    for u0 in u:
        d += (u0*u0)
    d = sqrt(d)
    return d
    
################################################
### gets the elements in a that are not in b ###
################################################
def setdiff(a,b):
    # INPUT
    #   - a: list with elements
    #   - b: list with elements
    # OUTPUT
    #   - aa: elements in a that are not in b
    b = set(b)
    return [aa for aa in a if aa not in b]

##########################################
#### gets all possible combinations   ####
#### for connection atoms in geometry ####
##########################################
def getbackbcombs():
    bbcombs = dict()
    bbcombs['one'] = [[1]]
    bbcombs['li'] = [[1],[2]]
    bbcombs['oct'] = [[1,2,3,4,5,6], # 6-dentate
                     [1,2,3,4,5],[1,2,3,4,6],[1,2,3,5,6],[1,2,4,5,6], # 5-dentate
                     [1,3,4,5,6],[2,3,4,5,6], # 5-dentate
                     [1,2,3,4],[2,5,4,6],[1,5,3,6], # 4-dentate
                     [1,2,3],[1,4,2],[1,4,3],[1,5,3],[1,6,3],[2,3,4], # 3-dentate
                     [2,5,4],[2,6,4],[5,4,6],[5,1,6],[5,2,6],[5,3,6], # 3-dentate
                     [1,2],[1,4],[1,5],[1,6],[2,3],[2,5], # 2-dentate
                     [2,6],[3,5],[3,6],[4,5],[4,6],[3,4], # 2-dentate
                     [1],[2],[3],[4],[5],[6]] # 1-dentate 
    bbcombs['pbp'] = [[1,2,3,4,5,6],[1,2,3,4,6], # 6/5-dentate
                      [1,2,3,5], # 4-dentate
                      [1,2,3],[1,2,4],[2,1,5],[3,1,6],[5,6,3],[2,6,5], # 3-dentate
                      [1,2],[2,3],[3,4],[4,5],[1,7],[2,6],[5,7],[3,6], # 2-dentate
                      [1],[2],[3],[4],[5],[6],[7]] # 1-dentate
    bbcombs['spy'] = [[1,2,3,4,5],[1,2,3,4],[1,2,3],[2,3,4],[3,4,1],[4,1,2],
                     [1,2],[1,4],[2,3],[3,4],[4,5],[2,5],[3,5],[1,5],[1],[2],[3],[4],[5]]
    bbcombs['sqp'] = [[1,4,2,3],[1,2,3],[2,3,4],[3,4,1],[4,1,2],[1,2],[1,4],[2,3],[3,4],
                      [1],[2],[3],[4]]
    bbcombs['tbp'] = [[1,2,3,4,5],[1,3,4,5],[3,2,4],[4,5,3],[5,1,3],[4,5],[5,3],[3,4],
                     [1,4],[1,5],[1,3],[2,4],[2,5],[2,3],[1],[2],[3],[4],[5]]
    bbcombs['thd'] = [[1,2,3,4],[3,2,4],[2,4,1],[4,1,3],[2,4],[4,3],[3,2],[1,3],[1,4],[2,4],[1],[2],[3],[4]]
    bbcombs['tpl'] = [[1,3,4],[1,2],[2,3],[1,3],[1],[2],[3]]
    bbcombs['tpr'] = [[1,2,3,4,5,6],[1,2,3,4,5],[1,2,5,4],[5,2,3,6],[1,4,6,3],[1,2,3],[3,6,5],
                     [2,3],[2,5],[5,6],[6,4],[4,1],[1],[2],[3],[4],[5],[6]]
    return bbcombs

##########################################
#### gets all possible combinations   ####
#### for connection atoms in geometry ####
####  in the case of forced order  #######
########   or unknown geometry   #########
##########################################
def getbackbcombsall(natoms):
    bbcombs = []
    nums = range(1,natoms+1)
    for i in range(1,natoms+1):
        bbcombs += list(itertools.combinations(nums,i))
    return bbcombs

#########################################
### gets a combination that satisfies ###
#### denticity and updates dictionary ###
#########################################
def getnupdateb(backbatoms,denticity):
    dlist = []
    batoms = []
    # find matching combination
    for b in backbatoms:
        if len(b)==denticity:
            batoms = b
            break
    # loop and find elements to delete
    for b in batoms:
        for i,bcomb in enumerate(backbatoms):
            if b in bcomb and i not in dlist:
                dlist.append(i)
    dlist.sort(reverse=True) # sort
    # delete used points
    for i in dlist:
        del backbatoms[i]
    if len(batoms) < 1:
        print 'No more connecting points available..'
    return batoms,backbatoms

##############################################
### gets connection atoms of smiles string ###
##############################################
def getsmilescat(args,indsmi):
    # INPUT
    #   - args: placeholder for input arguments
    #   - nosmiles: number of ligands defined via SMILES
    #   - indsmi: index of SMILES string ligand (like counter)
    # OUTPUT
    #   - tt: list of connection atoms
    tt= []  # initialize list of connection atoms
    if args.smicat and len(args.smicat)>indsmi: # get connection atom(s)
        tt = args.smicat[indsmi] # default value
    else:
        tt = [0] # default value 0 connection atom
    return tt

#######################################
### gets denticity of smiles string ###
#######################################
def getsmident(args,indsmi):
    # INPUT
    #   - args: placeholder for input arguments
    #   - nosmiles: number of ligands defined via SMILES
    #   - indsmi: index of SMILES string ligand (like counter)
    # OUTPUT
    #   - SMILES ligand denticity (int)
    ### check for denticity specification in input ###
    # if denticity is specified return this
    if args.smicat and len(args.smicat) > indsmi:
        return int(len(args.smicat[indsmi]))
    # otherwise return default
    else:
        return 1

##############################################
### modifies backbone according to pangles ###
##############################################
def modifybackbonep(backb, pangles):
    # INPUT
    #   - backb: list with points comprising the backbone
    #   - pangles: angles for distorting corresponding backbone points  (pairs of theta/phi)
    # OUTPUT
    #   - backb: list with modified points comprising the backbone
    for i,ll in enumerate(pangles):
        if ll:
            theta = pi*float(ll.split('/')[0])/180.0
            phi = pi*float(ll.split('/')[-1])/180.0
            backb[i+1] = PointTranslateSph(backb[0],backb[i+1],[distance(backb[0],backb[i+1]),theta,phi])
    return backb

##################################
### randomly distorts backbone ###
##################################
def distortbackbone(backb, distort):
    # INPUT
    #   - backb: list with points comprising the backbone
    #   - distort: % distortion of the backbone
    # OUTPUT
    #   - backb: list with modified points comprising the backbone
    for i in range(1,len(backb)):
            theta = random.uniform(0.0,0.01*int(distort)) # *0.5
            phi = random.uniform(0.0,0.01*int(distort)*0.5) # *0.5
            backb[i] = PointTranslateSph(backb[0],backb[i],[distance(backb[0],backb[i]),theta,phi])
    return backb
    
#######################
### reorder ligands ###
#######################
def smartreorderligs(args,ligs,dentl,licores):
    globs = globalvars()
    # INPUT
    #   - args: placeholder for input arguments
    #   - ligs: list of ligands
    #   - dents: ligand denticities
    # OUTPUT
    #   - indcs: reordering indices
    # check for forced order
    if not args.ligalign:
        indcs = range(0,len(ligs))
        return indcs
    lsizes = []
    for ligand in ligs:
        lig,emsg = lig_load(globs.installdir+'/',ligand,licores) # load ligand
        lsizes.append(len(lig.OBmol.atoms))
    # group by denticities
    dents = list(set(dentl))
    ligdentsidcs = [[] for a in dents]
    for i,dent in enumerate(dentl):
        ligdentsidcs[dents.index(dent)].append(i)
    # sort by highest denticity first
    ligdentsidcs = list(reversed(ligdentsidcs))
    indcs = []
    # within each group sort by size (smaller first)
    for ii,dd in enumerate(ligdentsidcs):
        locs = [lsizes[i] for i in dd]
        locind = [i[0] for i in sorted(enumerate(locs), key=lambda x:x[1])]
        for l in locind:
            indcs.append(ligdentsidcs[ii][l])
    return indcs
    
###########################################
### loads M-L bond length from database ###
###########################################
def getbondlength(args,metal,m3D,lig3D,matom,atom0,ligand,MLbonds):
    # INPUT
    #   - args: palceholder for input arguments
    #   - metal: name for metallic element
    #   - m3D: mol3D with main complex
    #   - lig3D: mol3D with ligand
    #   - matom: index of metal atom in m3D
    #   - atom0: index of connecting atom in lig3D
    #   - ligand: name of ligand
    #   - MLbonds: data from database
    # OUTPUT
    #   - bondl: bond length in A
    ### check for roman letters in oxstate
    romans={'I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6'}
    if args.oxstate: # if defined put oxstate in keys
        if args.oxstate in romans.keys():
            oxs = romans[args.oxstate]
        else:
            oxs = args.oxstate
    else:
        oxs = '-'
    # check for spin multiplicity
    spin = args.spin if args.spin else '-'
    key = []
    key.append((metal,oxs,spin,lig3D.getAtom(atom0).sym,ligand))
    key.append((metal,oxs,spin,lig3D.getAtom(atom0).sym,'-')) # disregard exact ligand
    key.append((metal,'-','-',lig3D.getAtom(atom0).sym,ligand)) # disregard oxstate/spin
    key.append((metal,'-','-',lig3D.getAtom(atom0).sym,'-')) # else just consider bonding atom
    found = False
    # search for data
    for kk in key:
        if (kk in MLbonds.keys()): # if exact key in dictionary
            bondl = float(MLbonds[kk])
            found = True
            break
    if not found: # last resort covalent radii
        bondl = m3D.getAtom(matom).rad + lig3D.getAtom(atom0).rad
    return bondl

###############################
### FORCE FIELD OPTIMIZATION ##
###############################
def ffopt(ff,mol,connected,constopt,frozenats,frozenangles,mlbonds):
    # INPUT
    #   - ff: force field to use, available MMFF94, UFF< Ghemical, GAFF
    #   - mol: mol3D to be ff optimized
    #   - connected: indices of connection atoms to metal
    #   - constopt: flag for constrained optimization
    # OUTPUT
    #   - mol: force field optimized mol3D
    ''' options
    B: ff before addition
    A: ff whole structure after addition (each step)
    BA: 0 and 1 together
    '''
    globs = globalvars()
    metals = range(21,31)+range(39,49)+range(72,81)
    ### check requested force field
    ffav = 'mmff94, uff, ghemical, gaff, mmff94s' # force fields
    if ff.lower() not in ffav:
        print 'Requested force field not available. Defaulting to MMFF94'
        ff = 'mmff94'
    # perform constrained ff optimization if requested after #
    if (constopt > 0):
        ### get metal
        midx = mol.findMetal()
        ### convert mol3D to OBmol via xyz file, because AFTER/END option have coordinates
        mol.writexyz(globs.homedir+'/tmp.xyz')
        mol.OBmol = mol.getOBmol(globs.homedir+'/tmp.xyz','xyzf')
        os.remove(globs.homedir+'/tmp.xyz')
        ### initialize constraints
        constr = openbabel.OBFFConstraints()
        ### openbabel indexing starts at 1 ### !!!
        # convert metals to carbons for FF
        indmtls = []
        mtlsnums = []
        for iiat,atom in enumerate(mol.OBmol.atoms):
            if atom.atomicnum in metals:
                indmtls.append(iiat)
                mtlsnums.append(atom.atomicnum)
                atom.OBAtom.SetAtomicNum(6)
        ### add distance constraints
        for ii,catom in enumerate(connected):
            if constopt==1 or frozenangles:
                constr.AddAtomConstraint(catom+1) # indexing babel
            else:
                constr.AddDistanceConstraint(midx+1,catom+1,mlbonds[ii]) # indexing babel
        for midxm in indmtls:
            ### freeze metal
            constr.AddAtomConstraint(midxm+1) # indexing babel
        ### freeze small ligands
        for cat in frozenats:
            constr.AddAtomConstraint(cat+1) # indexing babel
        ### set up forcefield
        forcefield = openbabel.OBForceField.FindForceField(ff)
        obmol = mol.OBmol.OBMol
        forcefield.Setup(obmol,constr)
        ### force field optimize structure
        if obmol.NumHvyAtoms() > 10:
            forcefield.ConjugateGradients(3000)
        else:
            forcefield.ConjugateGradients(2000)
        forcefield.GetCoordinates(obmol)
        en = forcefield.Energy()
        mol.OBmol = pybel.Molecule(obmol)
        # reset atomic number to metal
        for i,iiat in enumerate(indmtls):
            mol.OBmol.atoms[iiat].OBAtom.SetAtomicNum(mtlsnums[i])
        mol.convert2mol3D()
        del forcefield, constr, obmol
    else:
        ### initialize constraints
        constr = openbabel.OBFFConstraints()
        ### add atom constraints
        for catom in connected:
            constr.AddAtomConstraint(catom+1) # indexing babel
        ### set up forcefield
        forcefield = openbabel.OBForceField.FindForceField(ff)
        if len(connected) < 2:
            mol.OBmol.make3D('mmff94',1000) # add hydrogens and coordinates
        obmol = mol.OBmol.OBMol # convert to OBmol
        forcefield.Setup(obmol,constr)
        ### force field optimize structure
        if obmol.NumHvyAtoms() > 10:
            forcefield.ConjugateGradients(3000)
        else:
            forcefield.ConjugateGradients(2000)
        forcefield.GetCoordinates(obmol)
        en = forcefield.Energy()
        mol.OBmol = pybel.Molecule(obmol)
        mol.convert2mol3D()
        del forcefield, constr, obmol
    return mol,en
    
################################################
### FORCE FIELD OPTIMIZATION for custom cores ##
################################################
def ffoptd(ff,mol,connected,ccatoms,frozenats,nligats):
    # INPUT
    #   - ff: force field to use, available MMFF94, UFF< Ghemical, GAFF
    #   - mol: mol3D to be ff optimized
    #   - connected: indices of connection atoms to metal
    #   - constopt: flag for constrained optimization
    # OUTPUT
    #   - mol: force field optimized mol3D
    metals = range(21,31)+range(39,49)+range(72,81)
    ### convert mol3D to OBmol via xyz file, because AFTER/END option have coordinates
    mol.writexyz('tmp.xyz')
    mol.OBmol = mol.getOBmol('tmp.xyz','xyzf')
    os.remove('tmp.xyz')
    ### initialize constraints
    constr = openbabel.OBFFConstraints()
    ### openbabel indexing starts at 1 ### !!!
    # convert metals to carbons for FF
    indmtls = []
    mtlsnums = []
    for iiat,atom in enumerate(mol.OBmol.atoms):
        if atom.atomicnum in metals:
            indmtls.append(iiat)
            mtlsnums.append(atom.atomicnum)
            atom.OBAtom.SetAtomicNum(6)
    ### add distance constraints
    for ict,catom in enumerate(connected):
        dma = mol.getAtom(ccatoms[ict]).distance(mol.getAtom(catom))
        constr.AddDistanceConstraint(ccatoms[ict]+1,catom+1,dma) # indexing babel
    ### freeze core
    for ii in range(0,mol.natoms-nligats):
        constr.AddAtomConstraint(ii+1) # indexing babel
    ### freeze small ligands
    for cat in frozenats:
        constr.AddAtomConstraint(cat+1) # indexing babel
    ### set up forcefield
    forcefield = openbabel.OBForceField.FindForceField(ff)
    obmol = mol.OBmol.OBMol
    forcefield.Setup(obmol,constr)
    ### force field optimize structure
    if obmol.NumHvyAtoms() > 10:
        forcefield.ConjugateGradients(4000)
    else:
        forcefield.ConjugateGradients(2000)
    forcefield.GetCoordinates(obmol)
    en = forcefield.Energy()
    mol.OBmol = pybel.Molecule(obmol)
    # reset atomic number to metal
    for i,iiat in enumerate(indmtls):
        mol.OBmol.atoms[iiat].OBAtom.SetAtomicNum(mtlsnums[i])
    mol.convert2mol3D()
    del forcefield, constr, obmol
    return mol,en

#####################################################
####### Uses force field to estimate optimum ########
############ backbone positioning ###################
#################################################
def getconnection(core,cm,catom,toconnect):
    ff = 'UFF'
    metals = range(21,31)+range(39,49)+range(72,81)
    ### get hydrogens
    Hlist = core.getHs()
    ### add fake atoms for catoms
    ncore = core.natoms
    # add fake atom in local centermass axis
    coords = core.getAtom(catom).coords()
    dd = distance(coords,core.centermass())
    backbcoords = alignPtoaxis(coords,coords,vecdiff(coords,core.centermass()),1.5)
    bopt = backbcoords
    # manually find best positioning
    am = mol3D()
    am.addAtom(atom3D('C',backbcoords))
    for ii in range(0,toconnect-1):
        P = PointTranslateSph(coords,am.atoms[ii].coords(),[1.5,45,30])
        am.addAtom(atom3D('C',P))
    setopt = []
    mdist = -1
    for ii in range(0,toconnect):
        for itheta in range(0,360,3):
            for iphi in range(0,180,2):
                P = PointTranslateSph(coords,backbcoords,[1.5,itheta,iphi])
                am.atoms[ii].setcoords(P)
                dd = 0 
                for idx in range(0,toconnect):
                    dd += distance(cm,am.atoms[idx].coords())
                if (am.mindistmol() > 0.0):
                    d0 = dd+0.5*(log(core.mindist(am)*am.mindistmol()))
                if d0 > mdist:
                    mdist = d0
                    setopt = am.coordsvect()
    for ii in range(0,toconnect):
        core.addAtom(atom3D('C',setopt[ii]))
    ffoptc = False
    if ffoptc:
        ### convert mol3D to OBmol via xyz file
        core.writexyz('tmp.xyz')
        core.OBmol = core.getOBmol('tmp.xyz','xyzf')
        os.remove('tmp.xyz')
        ### openbabel indexing starts at 1 ### 
        # convert metals to carbons for FF
        [indmtls,mtlsnums] = [[],[]]
        for iiat,atom in enumerate(core.OBmol.atoms):
            if atom.atomicnum in metals:
                indmtls.append(iiat)
                mtlsnums.append(atom.atomicnum)
                core.OBmol.atoms[iiat].OBAtom.SetAtomicNum(6)
        ### initialize constraints
        constr = openbabel.OBFFConstraints()
        ### freeze molecule
        for atom in range(ncore):
            constr.AddAtomConstraint(atom+1) # indexing babel
        ### add distance constraints
        constr.AddDistanceConstraint(catom+1,ncore,1.5)
        ### ignore Hydrogens
        for ii in Hlist:
            constr.AddIgnore(ii+1)
        ### set up forcefield
        forcefield = openbabel.OBForceField.FindForceField(ff)
        obmol = core.OBmol.OBMol
        forcefield.Setup(obmol,constr)
        ### force field optimize structure
        forcefield.ConjugateGradients(500)
        forcefield.GetCoordinates(obmol)
        core.OBmol = pybel.Molecule(obmol)
        # reset atomic number to metal
        for i,iiat in enumerate(indmtls):
            core.OBmol.atoms[iiat].OBAtom.SetAtomicNum(mtlsnums[i])
        core.convert2mol3D()
        del forcefield, constr, obmol
    connPts = []
    for ii in range(0,toconnect):
        connPts.append(core.getAtom(ncore+ii).coords())
    return connPts
    
#################################################
####### functionalizes core with ligands ########
############## for metal complexes ##############
#################################################
def mcomplex(args,core,ligs,ligoc,installdir,licores,globs):
    # INPUT
    #   - args: placeholder for input arguments
    #   - core: mol3D structure with core
    #   - ligs: list of ligands
    #   - ligoc: list of ligand occupations
    #   - installdir: top installation directory
    #   - licores: dictionary with ligands
    #   - globs: class with global variables
    # OUTPUT
    #   - core3D: built complex
    #   - complex3D: list of all mol3D ligands and core
    #   - emsg: error messages
    if globs.debug:
        print '\nGenerating complex with ligands and occupations:',ligs,ligoc
    if args.gui:
        args.gui.iWtxt.setText('\n----------------------------------------------------------------------------------\n'+
                                      '\nGenerating complex with ligands: '+ ' '.join(ligs)+'\n'+args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    ### initialize variables ###
    emsg, complex3D = False, []
    # get available geometries
    coords,geomnames,geomshorts,geomgroups = getgeoms()
    coordbasef = geomgroups
    cclist = geomshorts # list of coordinations
    # get list of possible combinations for connectino atoms
    bbcombsdict = getbackbcombs()
    metal = core.getAtom(0).sym # metal symbol
    occs0 = []      # occurrences of each ligand
    cats0 = []      # connection atoms for ligands
    dentl = []      # denticity of ligands
    toccs = 0       # total occurrence count (number of ligands)
    octa = False    # flag for forced octahedral structures like porphyrines
    smilesligs = 0  # count how many smiles strings
    connected = []  # indices in core3D of ligand atoms connected to metal
    frozenats = []  # atoms to be frozen in optimization
    freezeangles = False # custom angles imposed
    MLoptbds = []   # list of bond lengths
    remCM = False   # remove dummy center of mass atom
    ### load bond data ###
    MLbonds = loaddata(installdir+'/Data/ML.dat')
    ### calculate occurrences, denticities etc for all ligands ###
    for i,ligname in enumerate(ligs):
        # if not in cores -> smiles/file
        if ligname not in licores.keys():
            if args.smicat and len(args.smicat) >= i and args.smicat[i]:
                if 'cm' in args.smicat[i]:
                    cats0.append(['c'])
                else:
                    cats0.append(args.smicat[i])
            else:
                cats0.append([1])
            dent_i = len(cats0[-1])
            smilesligs += 1
        else:
            cats0.append(False)
        # otherwise get denticity from ligands dictionary
            if 'cm' in licores[ligname][2]:
                dent_i = 1
            else:
                dent_i = int(len(licores[ligname][2]))
        # get occurrence for each ligand if specified (default 1)
        oc_i = int(ligoc[i]) if i < len(ligoc) else 1
        occs0.append(0)         # initialize occurrences list
        dentl.append(dent_i)    # append denticity to list
        # loop over occurrence of ligand i to check for max coordination
        for j in range(0,oc_i):
            occs0[i] += 1
            toccs += dent_i
    ### sort by descending denticity (needed for adjacent connection atoms) ###
    ligandsU,occsU,dentsU = ligs,occs0,dentl # save unordered lists
    indcs = smartreorderligs(args,ligs,dentl,licores)
    ligands = [ligs[i] for i in indcs]  # sort ligands list
    occs = [occs0[i] for i in indcs]    # sort occurrences list
    dents = [dentl[i] for i in indcs]   # sort denticities list
    tcats = [cats0[i] for i in indcs]# sort connections list
    # sort keepHs list ###
    keepHs = False
    if args.keepHs:
        keepHs = [k for k in args.keepHs]
        for j in range(len(args.keepHs),len(ligs)):
            keepHs.append(False)
        keepHs = [keepHs[i] for i in indcs] # sort keepHs list
    ### sort M-L bond list ###
    MLb = False
    if args.MLbonds:
        MLb = [k for k in args.MLbonds]
        for j in range(len(args.MLbonds),len(ligs)):
            MLb.append(False)
        MLb = [MLb[i] for i in indcs] # sort MLbonds list
    ### sort ligands custom angles ###
    pangles = False
    if args.pangles:
        pangles = []
        for j in range(len(args.pangles),len(ligs)):
            pangles.append(False)
        pangles = [args.pangles[i] for i in indcs] # sort custom langles list
    ### geometry information ###
    coord = toccs # complex coordination
    # check for coordination
    if args.coord and int(args.coord)!=coord:
        print "WARNING: Number of ligands doesn't agree with coordination/geometry. Will use geometry indicated by ligands."
        if args.gui:
            emsg = "Number of ligands doesn't agree with coordination/geometry. Will use geometry indicated by ligand frequency."
            qqb = mQDialogWarn('Warning',emsg)
            qqb.setParent(args.gui.wmain)
        if len(coordbasef) > coord -1 :
            geom = coordbasef[coord-1][0]
    elif args.coord:
        geom = coordbasef[int(args.coord)-1][0] # geometry specified by user coordination
    else:
        if len(coordbasef) > coord -1 :
            geom = coordbasef[coord-1][0] # total number of ligands define coordination
    # check if geometry is defined and overwrite
    if args.geometry and args.geometry in cclist:
        geom = args.geometry
    elif args.geometry:
        emsg = "Requested geometry not available."+"Defaulting to "+coordbasef[coord-1][0]
        if args.gui:
            qqb = mQDialogWarn('Warning',emsg)
            qqb.setParent(args.gui.wmain)
        print emsg
        print "Defaulting to "+coordbasef[coord-1][0]
    else:
        if len(coordbasef) <= coord-1:
            emsg = "WARNING: Coordination requested is not supproted. Defaulting to octahedral"
            print emsg
            if args.gui:
                qqb = mQDialogWarn('Warning',emsg)
                qqb.setParent(args.gui.wmain)
            geom = coordbasef[5][0] # force octahedrals
        else:
            geom = coordbasef[coord-1][0]
    ### load backbone and combinations ###
    # load backbone for coordination
    corexyz = loadcoord(installdir,geom)
    # get combinations possible for specified geometry
    if geom in bbcombsdict.keys() and not args.ligloc:
        backbatoms = bbcombsdict[geom]
    else:
        backbatoms = getbackbcombsall(len(corexyz)-1)
    # distort if requested
    if args.pangles:
        corexyz = modifybackbonep(corexyz,args.pangles) # point distortion
    if args.distort:
        corexyz = distortbackbone(corexyz,args.distort) # random distortion
    coord = len(corexyz)-1 # get coordination
    ### initialize molecules ###
    # create molecule and add metal and base
    m3D = mol3D() 
    m3D.addAtom(atom3D(metal,corexyz[0])) # add metal
    core3D = mol3D() # create backup
    core3D.addAtom(atom3D(metal,corexyz[0])) # add metal
    if args.calccharge:
        if args.oxstate:
            romans={'0':'0','I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6'}
            core3D.charge = int(romans[args.oxstate])
    mcoords = core3D.getAtom(0).coords() # metal coordinates in backbone
    ### initialize complex list of ligands/core
    auxm = mol3D()
    auxm.copymol3D(core3D)
    complex3D.append(auxm)
    # add terminal atoms in backbone given their coordinates
    for m in range(1,coord+1):
        m3D.addAtom(atom3D('X',corexyz[m])) ## add termination atoms
    #########################################################
    ####### Get connection points for all the ligands #######
    ########### smart alignment and forced order ############
    batslist = []
    if args.ligloc and args.ligalign:
        batslist0 = []
        for i,ligand in enumerate(ligandsU):
            for j in range(0,occsU[i]):
                # get correct atoms
                bats,backbatoms = getnupdateb(backbatoms,dentsU[i])
                batslist0.append(bats)
        # reorder according to smart reorder
        for i in indcs:
            offset = 0
            for ii in range(0,i):
                    offset += (occsU[ii]-1)
            for j in range(0,occsU[i]):
                batslist.append(batslist0[i+j+offset])# sort connections list
    else:
        for i,ligand in enumerate(ligands):
            for j in range(0,occs[i]):
                # get correct atoms
                bats,backbatoms = getnupdateb(backbatoms,dents[i])
                batslist.append(bats)
    #########################################################
    ###############################
    #### loop over ligands and ####
    ### begin functionalization ###
    ###############################
    # loop over ligands
    totlig = 0  # total number of ligands added
    ligsused = 0
    for i,ligand in enumerate(ligands):
        for j in range(0,occs[i]):
            denticity = dents[i]
            if not(ligand=='x' or ligand =='X') and (totlig-1+denticity < coord):
                # load ligand
                lig,emsg = lig_load(installdir,ligand,licores) # load ligand
                # check for smiles, force not removal of hydrogen
                allremH = True
                if ('+' in ligand or '-' in ligand):
                    allremH = False
                if emsg:
                    return False,emsg
                # if SMILES string
                if not lig.cat and tcats[i]:
                    if 'c' in tcats[i]:
                        lig.cat = [len(lig.OBmol.atoms)]
                    else:
                        lig.cat = tcats[i]
                # perform FF optimization if requested
                if args.ff and 'b' in args.ffoption:
                    if 'b' in lig.ffopt.lower():
                        lig,enl = ffopt(args.ff,lig,lig.cat,0,frozenats,freezeangles,MLoptbds)
                ###############################
                lig3D = lig # change name
                # convert to mol3D
                lig3D.convert2mol3D() # convert to mol3D
                if not keepHs or (len(keepHs) <= i or not keepHs[i]):
                    # remove one hydrogen
                    Hs = []
                    for cat in lig.cat:
                        Hs += lig3D.getHsbyIndex(cat)
                    if len(Hs) > 0 and allremH:
                        # check for cats indices
                        for ii,cat in enumerate(lig.cat):
                            if cat > Hs[0]:
                                lig.cat[ii] -= 1
                        lig3D.deleteatom(Hs[0])
                # check for CM connection atom
                if lig.cat[0] > lig3D.natoms-1 or 'cm' in lig.cat:
                    lig3D.addAtom(atom3D('C',lig3D.centermass()))
                    lig.cat = [lig3D.natoms-1]
                    remCM = True
                ### add atoms to connected atoms list
                catoms = lig.cat # connection atoms
                initatoms = core3D.natoms # initial number of atoms in core3D
                for at in catoms:
                    connected.append(initatoms+at)
                ### initialize variables
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0 # initialize variables
                ####################################################
                ##    attach ligand depending on the denticity    ##
                ## optimize geometry by minimizing steric effects ##
                ####################################################
                if (denticity == 1):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        emsg = 'Connecting all ligands is not possible. Check your input!'
                        if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    # connection atom in lig3D
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(lig3D.getAtom(atom0),m3D.getAtom(batoms[0]))
                    r1 = lig3D.getAtom(atom0).coords()
                    r2 = lig3D.centermass() # center of mass
                    if not r2:
                        emsg = 'Center of mass calculation for ligand failed. Check input.'
                        print emsg
                        if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    rrot = r1
                    theta,u = rotation_params(r0,r1,r2)
                    # for most ligands align center of mass of local environment
                    if (lig3D.natoms > 1):
                        lig3Db = mol3D()
                        lig3Db.copymol3D(lig3D)
                        ####################################
                        # center of mass of local environment (to avoid bad placement of bulky ligands)
                        auxmol = mol3D()
                        for at in lig3D.getBondedAtoms(atom0):
                            auxmol.addAtom(lig3D.getAtom(at))
                        r2 = auxmol.centermass() # overwrite global with local centermass
                        theta,u = rotation_params(r0,r1,r2)
                        ####################################
                        # rotate around axis and get both images
                        lig3D = rotate_around_axis(lig3D,rrot,u,theta)
                        lig3Db = rotate_around_axis(lig3Db,rrot,u,theta-180)
                        d2 = distance(mcoords,lig3D.centermass())
                        d1 = distance(mcoords,lig3Db.centermass())
                        lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    if lig3D.natoms > 2:
                        #####################################
                        # check for linear molecule
                        auxm = mol3D()
                        for at in lig3D.getBondedAtoms(atom0):
                            auxm.addAtom(lig3D.getAtom(at))
                        if auxm.natoms > 1:
                            r0 = lig3D.getAtom(atom0).coords()
                            r1 = auxm.getAtom(0).coords()
                            r2 = auxm.getAtom(1).coords()
                            if checkcolinear(r1,r0,r2):
                                # we will rotate so that angle is right
                                theta,urot = rotation_params(r1,mcoords,r2)
                                theta = vecangle(vecdiff(r0,mcoords),urot)
                                lig3D = rotate_around_axis(lig3D,r0,urot,theta)
                        #####################################
                        # check for symmetric molecule
                        if distance(lig3D.getAtom(atom0).coords(),lig3D.centersym()) < 8.0e-2:
                            atsc = lig3D.getBondedAtoms(atom0)
                            r0a = lig3D.getAtom(atom0).coords()
                            r1a = lig3D.getAtom(atsc[0]).coords()
                            r2a = lig3D.getAtom(atsc[1]).coords()
                            theta,u = rotation_params(r0a,r1a,r2a)
                            theta = vecangle(u,vecdiff(r0a,mcoords))
                            urot = cross(u,vecdiff(r0a,mcoords))
                            ####################################
                            # rotate around axis and get both images
                            lig3Db = mol3D()
                            lig3Db.copymol3D(lig3D)
                            lig3D = rotate_around_axis(lig3D,r0a,urot,theta)
                            lig3Db = rotate_around_axis(lig3Db,r0a,urot,-theta)
                            d2 = lig3D.mindist(core3D)
                            d1 = lig3Db.mindist(core3D)
                            lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                        # rotate around axis of symmetry and get best orientation
                        r1 = lig3D.getAtom(atom0).coords()
                        u = vecdiff(r1,mcoords)
                        dtheta = 2
                        optmax = -9999
                        totiters = 0
                        lig3Db = mol3D()
                        lig3Db.copymol3D(lig3D)
                        # check for minimum distance between atoms and center of mass distance
                        while totiters < 180:
                            lig3D = rotate_around_axis(lig3D,r1,u,dtheta)
                            d0 = lig3D.mindist(core3D) # try to maximize minimum atoms distance
                            d0cm = lig3D.distance(core3D) # try to maximize center of mass distance
                            iteropt = d0cm+10*log(d0) # optimization function
                            if (iteropt > optmax): # if better conformation, keep
                                lig3Db = mol3D()
                                lig3Db.copymol3D(lig3D)
                                optmax = iteropt
                            totiters += 1
                        lig3D = lig3Db
                    # get distance from bonds table or vdw radii
                    if MLb and MLb[i]:
                        if 'c' in MLb[i].lower():
                            bondl = m3D.getAtom(0).rad + lig3D.getAtom(atom0).rad
                        else:
                            bondl = float(MLb[i]) # check for custom
                    else:
                        bondl = getbondlength(args,metal,core3D,lig3D,0,atom0,ligand,MLbonds)
                    MLoptbds.append(bondl)
                    # get correct distance for center of mass
                    cmdist = bondl - distance(r1,mcoords)+distance(lig3D.centermass(),mcoords)
                    lig3D=setcmdistance(lig3D, mcoords, cmdist)
                elif (denticity == 2):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        if args.gui:
                            emsg = 'Connecting all ligands is not possible. Check your input!'
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    # connection atom
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(lig3D.getAtom(atom0),m3D.getAtom(batoms[0]))
                    r1 = lig3D.getAtom(atom0).coords()
                    # align center of mass to the middle
                    r21 = [a-b for a,b in zip(lig3D.getAtom(catoms[1]).coords(),r1)]
                    r21n = [a-b for a,b in zip(m3D.getAtom(batoms[1]).coords(),r1)]
                    if (norm(r21)*norm(r21n)) > 1e-8:
                        theta = 180*arccos(dot(r21,r21n)/(norm(r21)*norm(r21n)))/pi
                    else:
                        theta = 0.0
                    u = cross(r21,r21n)
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,r1,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,r1,u,theta-180)
                    d1 = distance(lig3D.getAtom(catoms[1]).coords(),m3D.getAtom(batoms[1]).coords())
                    d2 = distance(lig3Db.getAtom(catoms[1]).coords(),m3D.getAtom(batoms[1]).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # flip if overlap
                    r0l = lig3D.getAtom(catoms[0]).coords()
                    r1l = lig3D.getAtom(catoms[1]).coords()
                    md = min(distance(r0l,mcoords),distance(r1l,mcoords))
                    if lig3D.mindist(core3D) < md:
                        lig3D = rotate_around_axis(lig3D,r0l,vecdiff(r1l,r0l),180.0)
                    # correct plane
                    r0b = m3D.getAtom(batoms[0]).coords()
                    r1b = m3D.getAtom(batoms[1]).coords()
                    r0l = lig3D.getAtom(catoms[0]).coords()
                    r1l = lig3D.getAtom(catoms[1]).coords()
                    rm = lig3D.centermass()
                    urot = vecdiff(r1l,r0l)
                    theta,ub = rotation_params(mcoords,r0b,r1b)
                    theta,ul = rotation_params(rm,r0l,r1l)
                    if (norm(ub)*norm(ul)) > 1e-8:
                        theta = 180*arccos(dot(ub,ul)/(norm(ub)*norm(ul)))/pi-180.0
                    else:
                        theta = 0.0
                    # rotate around axis 
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    lig3D = rotate_around_axis(lig3D,r1,urot,theta)
                    lig3Db = rotate_around_axis(lig3Db,r1,urot,-theta)
                    # select best
                    rm0,rm1 = lig3D.centermass(),lig3Db.centermass()
                    theta,ul0 = rotation_params(rm0,r0l,r1l)
                    theta,ul1 = rotation_params(rm1,r0l,r1l)
                    th0 = 180*arccos(dot(ub,ul0)/(norm(ub)*norm(ul0)))/pi
                    th0 = min(abs(th0),abs(180-th0))
                    th1 = 180*arccos(dot(ub,ul1)/(norm(ub)*norm(ul1)))/pi
                    th1 = min(abs(th1),abs(180-th1))
                    lig3D = lig3D if th0 < th1 else lig3Db
                    # get distance from bonds table or vdw radii
                    if MLb and MLb[i]:
                        if 'c' in MLb[i].lower():
                            bondl = m3D.getAtom(0).rad + lig3D.getAtom(atom0).rad
                        else:
                            bondl = float(MLb[i]) # check for custom
                    else:
                        bondl = getbondlength(args,metal,core3D,lig3D,0,atom0,ligand,MLbonds)
                    MLoptbds.append(bondl)
                    MLoptbds.append(bondl)
                    lig3D = setPdistance(lig3D, r1, r0, bondl)
                    # fix ML bond length and distort angle if needed
                    rtarget = getPointu(mcoords, bondl, vecdiff(r1b,mcoords)) # get second point target
                    dr = vecdiff(rtarget,lig3D.getAtom(catoms[1]).coords())
                    # distort ligand in nsteps steps
                    nsteps = 15 
                    ddr = [di/nsteps for di in dr]
                    ens =[]
                    cutoff = 5.0 # kcal/mol
                    for ii in range(0,nsteps):
                        lig3D,enl = ffopt('mmff94',lig3D,[],1,[catoms[0],catoms[1]],False,[])
                        ens.append(enl)
                        lig3D.getAtom(catoms[1]).translate(ddr)
                        # check fo cutoff
                        if ens[-1] - ens[0] > 5.0:
                            # fix ML bond length get optimum guess
                            r0,r1 = lig3D.getAtomCoords(catoms[0]),lig3D.getAtomCoords(catoms[1])
                            r01 = distance(r0,r1)
                            theta1 = 180*arccos(0.5*r01/bondl)/pi
                            theta2 = vecangle(vecdiff(r1,r0),vecdiff(mcoords,r0))
                            dtheta = theta2-theta1
                            theta,urot = rotation_params(mcoords,r0,r1)
                            lig3D = rotate_around_axis(lig3D,r0,urot,-dtheta) # rotate so that it matches bond
                            break
                    # freeze local geometry
                    lats = lig3D.getBondedAtoms(catoms[0])+lig3D.getBondedAtoms(catoms[1])
                    for lat in list(set(lats)):
                        frozenats.append(lat+core3D.natoms)
                elif (denticity == 3):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        if args.gui:
                            emsg = 'Connecting all ligands is not possible. Check your input!'
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    # connection atom
                    atom0 = catoms[1]
                    ### align molecule according to connection atom and shadow atom ###
                    lig3D.alignmol(lig3D.getAtom(atom0),m3D.getAtom(batoms[1]))
                    # align with correct plane
                    rl0,rl1,rl2 = lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),lig3D.getAtom(catoms[2]).coords()
                    rc0,rc1,rc2 = m3D.getAtom(batoms[0]).coords(),m3D.getAtom(batoms[1]).coords(),m3D.getAtom(batoms[2]).coords()
                    theta0,ul = rotation_params(rl0,rl1,rl2)
                    theta1,uc = rotation_params(rc0,rc1,rc2)
                    urot = vecdiff(rl1,mcoords)
                    theta = vecangle(ul,uc)
                    ### rotate around primary axis ###
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    lig3D = rotate_around_axis(lig3D,rl1,urot,theta)
                    lig3Db = rotate_around_axis(lig3Db,rl1,urot,180-theta)
                    rl0,rl1,rl2 = lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),lig3D.getAtom(catoms[2]).coords()
                    rl0b,rl1b,rl2b = lig3Db.getAtom(catoms[0]).coords(),lig3Db.getAtom(catoms[1]).coords(),lig3Db.getAtom(catoms[2]).coords()
                    rc0,rc1,rc2 = m3D.getAtom(batoms[0]).coords(),m3D.getAtom(batoms[1]).coords(),m3D.getAtom(batoms[2]).coords()
                    theta,ul = rotation_params(rl0,rl1,rl2)
                    theta,ulb = rotation_params(rl0b,rl1b,rl2b)
                    theta,uc = rotation_params(rc0,rc1,rc2)
                    d1 = norm(cross(ul,uc))
                    d2 = norm(cross(ulb,uc))
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    ### rotate around secondary axis ###
                    auxm = mol3D()
                    auxm.addAtom(lig3D.getAtom(catoms[0]))
                    auxm.addAtom(lig3D.getAtom(catoms[2]))
                    theta,urot0 = rotation_params(core3D.getAtom(0).coords(),lig3D.getAtom(atom0).coords(),auxm.centermass())
                    theta0,urot = rotation_params(lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),lig3D.getAtom(catoms[2]).coords())
                    # change angle if > 90
                    if theta > 90:
                        theta -= 180
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    lig3D = rotate_around_axis(lig3D,lig3D.getAtom(atom0).coords(),urot,theta)
                    lig3Db = rotate_around_axis(lig3Db,lig3D.getAtom(atom0).coords(),urot,180-theta)
                    d1 = distance(lig3D.getAtom(catoms[0]).coords(),m3D.getAtom(batoms[0]).coords())
                    d2 = distance(lig3Db.getAtom(catoms[0]).coords(),m3D.getAtom(batoms[0]).coords())
                    lig3D = lig3D if (d1 < d2) else lig3Db
                    # correct if not symmetric
                    theta0,urotaux = rotation_params(lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),core3D.getAtom(0).coords())
                    theta1,urotaux = rotation_params(lig3D.getAtom(catoms[2]).coords(),lig3D.getAtom(catoms[1]).coords(),core3D.getAtom(0).coords())
                    dtheta = 0.5*(theta1-theta0)
                    if abs(dtheta) > 0.5:
                        lig3D = rotate_around_axis(lig3D,lig3D.getAtom(atom0).coords(),urot,dtheta)
                    # flip to align 3rd atom if wrong
                    urot = vecdiff(lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords())
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    lig3Db = rotate_around_axis(lig3Db,rc1,urot,180)
                    d1 = distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords())
                    d2 = distance(lig3Db.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    # if overlap flip
                    dm0 = distance(lig3D.getAtom(catoms[0]).coords(),m3D.getAtom(0).coords())
                    dm1 = distance(lig3D.getAtom(catoms[1]).coords(),m3D.getAtom(0).coords())
                    dm2 = distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(0).coords())
                    mind = min([dm0,dm1,dm2])
                    for iiat,atom in enumerate(lig3D.atoms):
                        if iiat not in catoms and distance(atom.coords(),m3D.getAtom(0).coords()) < min([dm0,dm1,dm2]):
                            lig3D = rotate_around_axis(lig3D,rc1,uc,180)
                            break
                    # get distance from bonds table or vdw radii
                    if MLb and MLb[i]:
                        if 'c' in MLb[i].lower():
                            bondl = m3D.getAtom(0).rad + lig3D.getAtom(atom0).rad
                        else:
                            bondl = float(MLb[i]) # check for custom
                    else:
                        bondl = getbondlength(args,metal,core3D,lig3D,0,atom0,ligand,MLbonds)
                    for iib in range(0,3):
                        MLoptbds.append(bondl)
                    # set correct distance
                    setPdistance(lig3D, lig3D.getAtom(atom0).coords(), m3D.getAtom(0).coords(), bondl)
                elif (denticity == 4):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        if args.gui:
                            emsg = 'Connecting all ligands is not possible. Check your input!'
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    # connection atom
                    atom0 = catoms[0]
                    # align molecule according to symmetry center
                    auxmol = mol3D()
                    for iiax in range(0,4):
                        auxmol.addAtom(lig3D.getAtom(catoms[iiax]))
                    lig3D.alignmol(atom3D('C',auxmol.centermass()),m3D.getAtom(0))
                    # align plane
                    r0c = m3D.getAtom(batoms[0]).coords()
                    r1c = m3D.getAtom(batoms[1]).coords()
                    r2c = m3D.getAtom(batoms[2]).coords()
                    r0l = lig3D.getAtom(catoms[0]).coords()
                    r1l = lig3D.getAtom(catoms[1]).coords()
                    r2l = lig3D.getAtom(catoms[2]).coords()
                    theta,uc = rotation_params(r0c,r1c,r2c) # normal vector to backbone plane
                    theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    theta = 180*arccos(dot(uc,ul)/(norm(uc)*norm(ul)))/pi
                    u = cross(uc,ul)
                    # rotate around axis to match planes
                    theta = 180-theta if theta > 90 else theta
                    lig3D = rotate_around_axis(lig3D,r0l,u,theta)
                    # rotate around secondary axis to match atoms
                    r0l = lig3D.getAtom(catoms[0]).coords()
                    r1l = lig3D.getAtom(catoms[1]).coords()
                    r2l = lig3D.getAtom(catoms[2]).coords()
                    theta0,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    rm = lig3D.centermass()
                    r1 = vecdiff(r0l,mcoords)
                    r2 = vecdiff(r0c,mcoords)
                    theta = 180*arccos(dot(r1,r2)/(norm(r1)*norm(r2)))/pi
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    # rotate around axis and get both images
                    lig3D = rotate_around_axis(lig3D,mcoords,ul,theta)
                    # get distance from bonds table or vdw radii
                    if MLb and MLb[i]:
                        if 'c' in MLb[i].lower():
                            bondl = m3D.getAtom(0).rad + lig3D.getAtom(atom0).rad
                        else:
                            bondl = float(MLb[i]) # check for custom
                    else:
                        bondl = getbondlength(args,metal,core3D,lig3D,0,atom0,ligand,MLbonds)
                    for iib in range(0,4):
                        MLoptbds.append(bondl)
                elif (denticity == 5):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        emsg = 'Connecting all ligands is not possible. Check your input!'
                        break
                    # get center of mass 
                    ligc = mol3D()
                    for i in range(0,4): #5 is the non-planar atom
                        ligc.addAtom(lig3D.getAtom(catoms[i]))
                    # translate ligand to the middle of octahedral
                    lig3D.translate(vecdiff(mcoords,ligc.centermass()))
                    # get plane
                    r0c = m3D.getAtom(batoms[0]).coords()
                    r2c = m3D.getAtom(batoms[1]).coords()
                    r1c = mcoords
                    r0l = lig3D.getAtom(catoms[0]).coords()
                    r2l = lig3D.getAtom(catoms[1]).coords()
                    r1l = mcoords
                    theta,uc = rotation_params(r0c,r1c,r2c) # normal vector to backbone plane
                    theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    theta = vecangle(uc,ul)
                    u = cross(uc,ul)
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    # rotate around axis to match planes
                    lig3D = rotate_around_axis(lig3D,mcoords,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,mcoords,u,180+theta)
                    d1 = distance(lig3D.getAtom(catoms[4]).coords(),m3D.getAtom(batoms[-1]).coords())
                    d2 = distance(lig3Db.getAtom(catoms[4]).coords(),m3D.getAtom(batoms[-1]).coords())
                    lig3D = lig3D if (d2 < d1)  else lig3Db # pick best one
                    # rotate around center axis to match backbone atoms
                    r0l = vecdiff(lig3D.getAtom(catoms[0]).coords(),mcoords)
                    r1l = vecdiff(m3D.getAtom(totlig+1).coords(),mcoords)
                    u = cross(r0l,r1l)
                    theta = 180*arccos(dot(r0l,r1l)/(norm(r0l)*norm(r1l)))/pi
                    lig3Db = mol3D()
                    lig3Db.copymol3D(lig3D)
                    lig3D = rotate_around_axis(lig3D,mcoords,u,theta)
                    lig3Db = rotate_around_axis(lig3Db,mcoords,u,theta-90)
                    d1 = distance(lig3D.getAtom(catoms[0]).coords(),m3D.getAtom(batoms[0]).coords())
                    d2 = distance(lig3Db.getAtom(catoms[0]).coords(),m3D.getAtom(batoms[0]).coords())
                    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    bondl = getbondlength(args,metal,core3D,lig3D,0,catoms[0],ligand,MLbonds)
                    # flip if necessary
                    if len(batslist) > ligsused:
                        nextatbats = batslist[ligsused]
                    auxm = mol3D()
                    if len(nextatbats) > 0:
                        for at in nextatbats:
                            auxm.addAtom(m3D.getAtom(at))
                        if lig3D.overlapcheck(auxm,True): # if overlap flip
                            urot = vecdiff(m3D.getAtomCoords(batoms[1]),m3D.getAtomCoords(batoms[0]))
                            lig3D = rotate_around_axis(lig3D,mcoords,urot,180)
                    for iib in range(0,5):
                        MLoptbds.append(bondl)
                elif (denticity == 6):
                    # connection atoms in backbone
                    batoms = batslist[ligsused]
                    if len(batoms) < 1 :
                        if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        emsg = 'Connecting all ligands is not possible. Check your input!'
                        break
                    # get center of mass 
                    ligc = mol3D()
                    for i in range(0,6):
                        ligc.addAtom(lig3D.getAtom(catoms[i]))
                    # translate metal to the middle of octahedral
                    core3D.translate(vecdiff(ligc.centermass(),mcoords))
                    bondl = getbondlength(args,metal,core3D,lig3D,0,catoms[0],ligand,MLbonds)
                    for iib in range(0,6):
                        MLoptbds.append(bondl)
                auxm = mol3D()
                auxm.copymol3D(lig3D)
                complex3D.append(auxm)
                if 'a' not in lig.ffopt.lower():
                    for latdix in range(0,lig3D.natoms):
                        frozenats.append(latdix+core3D.natoms)
                # combine molecules
                core3D = core3D.combine(lig3D)
                if args.calccharge:
                    core3D.charge += lig3D.charge
                # remove dummy cm atom if requested
                if remCM:
                    core3D.deleteatom(core3D.natoms-1)
                # perform FF optimization if requested
                if args.ff and 'a' in args.ffoption:
                    core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds)
            totlig += denticity
            ligsused += 1
    # perform FF optimization if requested
    if args.ff and 'a' in args.ffoption:
        core3D,enc = ffopt(args.ff,core3D,connected,2,frozenats,freezeangles,MLoptbds)
    ###############################
    return core3D,complex3D,emsg

#################################################
####### functionalizes core with ligands ########
############## for metal complexes ##############
#################################################
def customcore(args,core,ligs,ligoc,installdir,licores,globs):
    # INPUT
    #   - args: placeholder for input arguments
    #   - core: mol3D structure with core
    #   - ligs: list of ligands
    #   - ligoc: list of ligand occupations
    #   - installdir: top installation directory
    #   - licores: dictionary with ligands
    #   - globs: class with global variables
    # OUTPUT
    #   - core3D: built complex
    #   - complex3D: list of all mol3D ligands and core
    #   - emsg: error messages
    if globs.debug:
        print '\nGenerating complex with ligands and occupations:',ligs,ligoc
    if args.gui:
        args.gui.iWtxt.setText('\nGenerating complex with core:'+args.core+' and ligands: '+ ' '.join(ligs)+'\n'+args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    ### initialize variables ###
    emsg, complex3D = False, []
    occs0 = []      # occurrences of each ligand
    toccs = 0       # total occurrence count (number of ligands)
    catsmi = []     # SMILES ligands connection atoms
    smilesligs = 0  # count how many smiles strings
    cats0 = []
    dentl = []      # denticity of ligands
    connected = []  # indices in core3D of ligand atoms connected to metal
    frozenats = []  # list of frozen atoms for optimization
    ### load bond data ###
    MLbonds = loaddata(installdir+'/Data/ML.dat')
    ### calculate occurrences, denticities etc for all ligands ###
    for i,ligname in enumerate(ligs):
        # if not in cores -> smiles/file
        if ligname not in licores.keys():
            if args.smicat and len(args.smicat) >= i and args.smicat[i]:
                cats0.append(args.smicat[i])
            else:
                cats0.append([1])
            dent_i = len(cats0[-1])
            smilesligs += 1
        else:
            cats0.append(False)
        # otherwise get denticity from ligands dictionary
            dent_i = int(len(licores[ligname][2]))
        # get occurrence for each ligand if specified (default 1)
        oc_i = int(ligoc[i]) if i < len(ligoc) else 1
        occs0.append(0)         # initialize occurrences list
        dentl.append(dent_i)    # append denticity to list
        # loop over occurrence of ligand i to check for max coordination
        for j in range(0,oc_i):
            if (toccs+dent_i <= 7):
                occs0[i] += 1
            toccs += dent_i
    # remove ligands with denticity > 1
    todel = []
    for ii,ddent in enumerate(dentl):
        if ddent > 1:
            todel.append(ii)
    for ii in sorted(todel,reverse=True):
        del dentl[ii]
        del ligands[ii]
        del occs[ii]
    ### sort by descending denticity (needed for adjacent connection atoms) ###
    indcs = smartreorderligs(args,ligs,dentl,licores)
    ligands = [ligs[i] for i in indcs]  # sort ligands list
    occs = [occs0[i] for i in indcs]    # sort occurrences list
    tcats = [cats0[i] for i in indcs]# sort issmiles list
    # sort keepHs list ###
    keepHs = False
    if args.keepHs:
        keepHs = [k for k in args.keepHs]
        for j in range(len(args.keepHs),len(ligs)):
            keepHs.append(False)
        keepHs = [keepHs[i] for i in indcs] # sort keepHs list
    ### sort M-L bond list ###
    MLb = False
    if args.MLbonds:
        MLb = [k for k in args.MLbonds]
        for j in range(len(args.MLbonds),len(ligs)):
            MLb.append(False)
        MLb = [MLb[i] for i in indcs] # sort MLbonds list
    if not args.ccatoms:
        emsg = 'Connection atoms for custom core not specified. Defaulting to 1!\n'
        print emsg
        if args.gui:
            qqb = mQDialogWarn('Warning',emsg)
            qqb.setParent(args.gui.wmain)
    ccatoms = args.ccatoms if args.ccatoms else [0]
    core3D = mol3D()
    core3D.copymol3D(core)
    cmcore = core3D.centermass()
    if args.calccharge:
        if args.oxstate:
            romans={'0':'0','I':'1','II':'2','III':'3','IV':'4','V':'5','VI':'6'}
            core3D.charge = int(romans[args.oxstate])
    # remove one hydrogen for each functionalization
    Hs = []
    if not args.replig:
        for ccat in ccatoms:
            Hs += core3D.getHsbyAtom(core3D.getAtom(ccat))
    # remove hydrogens and shift ccatoms
    if len(Hs) > 0:
        for H in sorted(Hs,reverse=True):
            core3D.deleteatom(H)
            # fix indexing
            for ii,cat in enumerate(ccatoms):
                if cat > H:
                    ccatoms[ii] -= 1
    ###############################
    #### loop over ligands and ####
    ### begin functionalization ###
    ###############################
    # flags of connection points already calculated
    setccatoms = list(set(ccatoms)) # set of connection points
    conflags = [False for ii in enumerate(setccatoms)]
    confcount = 0 
    # loop over ligands
    totlig = 0  # total number of ligands added
    for i,ligand in enumerate(ligands):
        if len(ccatoms) < i:
            ccatoms.append(0)
        for j in range(0,occs[i]):
            if not(ligand=='x' or ligand =='X'):
                if totlig >= len(ccatoms):
                    emsg = 'Number of ligands greater than connection points. Please specify enough connection atoms in custom core.\n'
                    print emsg
                    if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                    return False,emsg
                core = mol3D()
                core.copymol3D(core3D)
                if not args.replig:
                    cidxconn = setccatoms.index(ccatoms[totlig]) # get current index in set of connection atoms
                    if not conflags[cidxconn]:
                        totconn = ccatoms.count(ccatoms[totlig]) # total connection points
                        alconn = ccatoms[:totlig].count(ccatoms[totlig]) # already connected
                        cpoints = getconnection(core,cmcore,ccatoms[totlig],totconn)
                        conflags[cidxconn] = True
                        confcount = 0 
                    else:
                        confcount += 1
                    cpoint = cpoints[confcount]
                    mcoords = core3D.getAtom(ccatoms[totlig]).coords() # metal coordinates in backbone
                    # connection atom save
                    conatom3D = atom3D(core3D.getAtom(ccatoms[totlig]).sym,core3D.getAtom(ccatoms[totlig]).coords())
                else:
                    cpoint = core3D.getAtom(ccatoms[totlig]).coords()
                    conatoms = core3D.getBondedAtoms(ccatoms[totlig])
                    # find smaller ligand to remove
                    minmol = 10000
                    mindelats = []
                    atclose = 0
                    # loop over different connected atoms
                    for cat in conatoms:
                        # find submolecule
                        delatoms = core3D.findsubMol(ccatoms[totlig],cat) 
                        if len(delatoms) < minmol: # check for smallest
                            mindelats = delatoms
                            minmol = len(delatoms) # size
                            atclose = cat # connection atom
                        # if same atoms in ligand get shortest distance
                        elif len(delatoms)==minmol:
                            d0 = core3D.getAtom(ccatoms[totlig]).distance(core3D.getAtom(cat))
                            d1 = core3D.getAtom(ccatoms[totlig]).distance(core3D.getAtom(mindelats[0]))
                            if d0 < d1:
                                mindelats = delatoms
                                atclose = cat
                    mcoords = core3D.getAtom(atclose).coords() # connection coordinates in backbone
                    # connection atom save
                    conatom3D = atom3D(core3D.getAtom(atclose).sym,core3D.getAtom(atclose).coords())
                    delatoms = mindelats
                    # find shifting if needed
                    if len(ccatoms) > totlig+1:
                        for cccat in range(totlig+1,len(ccatoms)):
                            lshift = len([a for a in delatoms if a < ccatoms[cccat]])
                            ccatoms[cccat] -= lshift
                    core3D.deleteatoms(delatoms)
                # check for smiles, force not removal of hydrogen
                allremH = True
                if ('+' in ligand or '-' in ligand):
                    allremH = False
                # load ligand
                lig,emsg = lig_load(installdir,ligand,licores) # load ligand
                if emsg:
                    return False,emsg
                # if SMILES string
                if not lig.cat and tcats[i]:
                    lig.cat = tcats[i]
                # perform FF optimization if requested
                if args.ff and 'b' in args.ffoption:
                    if 'B' in lig.ffopt:
                        lig,enl = ffopt(args.ff,lig,lig.cat,0,frozenats,False,False)
                ###############################
                lig3D = lig # change name
                # convert to mol3D
                lig3D.convert2mol3D() # convert to mol3D
                if not keepHs or (len(keepHs) <= i or not keepHs[i]):
                    # remove one hydrogen
                    Hs = lig3D.getHsbyIndex(lig.cat[0])
                    if len(Hs) > 0 and allremH:
                        lig3D.deleteatom(Hs[0])
                ### add atoms to connected atoms list
                catoms = lig.cat # connection atoms
                initatoms = core3D.natoms # initial number of atoms in core3D
                ### initialize variables
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0 # initialize variables
                ####################################################
                ##    attach ligand depending on the denticity    ##
                ## optimize geometry by minimizing steric effects ##
                ####################################################
                denticity = 1
                if (denticity == 1):
                    # connection atom in lig3D
                    atom0 = catoms[0]
                    # align molecule according to connection atom and shadow atom
                    lig3D.alignmol(lig3D.getAtom(atom0),atom3D('C',cpoint))
                    r1 = lig3D.getAtom(atom0).coords()
                    r2 = lig3D.centermass() # center of mass
                    if not r2:
                        emsg = 'Center of mass calculation for ligand failed. Check input.'
                        print emsg
                        if args.gui:
                            qqb = mQDialogWarn('Warning',emsg)
                            qqb.setParent(args.gui.wmain)
                        break
                    rrot = r1
                    theta,u = rotation_params(r0,r1,r2)
                    # for most ligands align center of mass of local environment
                    if (lig3D.natoms > 1):
                        lig3Db = mol3D()
                        lig3Db.copymol3D(lig3D)
                        ####################################
                        # center of mass of local environment (to avoid bad placement of bulky ligands)
                        auxmol = mol3D()
                        for at in lig3D.getBondedAtoms(atom0):
                            auxmol.addAtom(lig3D.getAtom(at))
                        r2 = auxmol.centermass() # overwrite global with local centermass
                        theta,u = rotation_params(r0,r1,r2)
                        ####################################
                        # rotate around axis and get both images
                        lig3D = rotate_around_axis(lig3D,rrot,u,theta)
                        lig3Db = rotate_around_axis(lig3Db,rrot,u,theta-180)
                        d2 = distance(mcoords,lig3D.centermass())
                        d1 = distance(mcoords,lig3Db.centermass())
                        lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                    if lig3D.natoms > 2:
                        #####################################
                        # check for linear molecule
                        auxm = mol3D()
                        for at in lig3D.getBondedAtoms(atom0):
                            auxm.addAtom(lig3D.getAtom(at))
                        if auxm.natoms > 1:
                            r0 = lig3D.getAtom(atom0).coords()
                            r1 = auxm.getAtom(0).coords()
                            r2 = auxm.getAtom(1).coords()
                            if checkcolinear(r1,r0,r2):
                                # we will rotate so that angle is right
                                theta,urot = rotation_params(r1,mcoords,r2)
                                theta = vecangle(vecdiff(r0,mcoords),urot)
                                lig3D = rotate_around_axis(lig3D,r0,urot,theta)
                        #####################################
                        # check for symmetric molecule
                        if distance(lig3D.getAtom(atom0).coords(),lig3D.centersym()) < 8.0e-2:
                            atsc = lig3D.getBondedAtoms(atom0)
                            r0a = lig3D.getAtom(atom0).coords()
                            r1a = lig3D.getAtom(atsc[0]).coords()
                            r2a = lig3D.getAtom(atsc[1]).coords()
                            theta,u = rotation_params(r0a,r1a,r2a)
                            theta = vecangle(u,vecdiff(r0a,mcoords))
                            urot = cross(u,vecdiff(r0a,mcoords))
                            ####################################
                            # rotate around axis and get both images
                            lig3Db = mol3D()
                            lig3Db.copymol3D(lig3D)
                            lig3D = rotate_around_axis(lig3D,r0a,urot,theta)
                            lig3Db = rotate_around_axis(lig3Db,r0a,urot,-theta)
                            d2 = lig3D.mindist(core3D)
                            d1 = lig3Db.mindist(core3D)
                            lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
                        # rotate around axis of symmetry and get best orientation
                        r1 = lig3D.getAtom(atom0).coords()
                        u = vecdiff(r1,mcoords)
                        dtheta = 2
                        optmax = -9999
                        totiters = 0
                        lig3Db = mol3D()
                        lig3Db.copymol3D(lig3D)
                        # check for minimum distance between atoms and center of mass distance
                        while totiters < 180:
                            lig3D = rotate_around_axis(lig3D,r1,u,dtheta)
                            d0 = lig3D.mindist(core3D) # try to maximize minimum atoms distance
                            d0cm = lig3D.distance(core3D) # try to maximize center of mass distance
                            iteropt = d0cm+10*log(d0) # optimization function
                            if (iteropt > optmax): # if better conformation, keep
                                lig3Db = mol3D()
                                lig3Db.copymol3D(lig3D)
                                optmax = iteropt
                            totiters += 1
                        lig3D = lig3Db
                    # get distance from bonds table or vdw radii
                    if MLb and MLb[i]:
                        if 'c' in MLb[i].lower():
                            bondl = conatom3D.rad + lig3D.getAtom(atom0).rad
                        else:
                            bondl = float(MLb[i]) # check for custom
                    else:
                        mm3D = mol3D()
                        mm3D.addAtom(conatom3D)
                        bondl = getbondlength(args,conatom3D.sym,mm3D,lig3D,0,atom0,ligand,MLbonds)
                    # get correct distance for center of mass
                    u = vecdiff(cpoint,mcoords)
                    lig3D = aligntoaxis2(lig3D, cpoint, mcoords, u, bondl)
                    connected.append(core3D.natoms+atom0)
                    # list of frozen atoms (small ligands)
                    if 'A' not in lig.ffopt:
                        for latdix in range(0,lig3D.natoms):
                            frozenats.append(latdix+core3D.natoms)
                    # combine molecules
                    core3D = core3D.combine(lig3D)
                else:
                    emsg = 'Multidentate ligands not supported for custom cores. Skipping.\n' 
                    print emsg
                if args.calccharge:
                    core3D.charge += lig3D.charge
                nligats = lig3D.natoms
                # perform FF optimization if requested
                if args.ff and 'a' in args.ffoption:
                    core3D,enc = ffoptd(args.ff,core3D,connected,ccatoms,frozenats,nligats)
            totlig += 1
    # perform FF optimization if requested
    if args.ff and 'a' in args.ffoption:
        core3D,enc = ffoptd(args.ff,core3D,connected,ccatoms,frozenats,nligats)
    return core3D,emsg

##########################################
### main structure generation function ###
##########################################
def structgen(installdir,args,rootdir,ligands,ligoc,globs):
    # INPUT
    #   - installdir: top installation directory
    #   - args: placeholder for input arguments
    #   - rootdir: directory of current run
    #   - ligands: list of ligands
    #   - ligoc: list of ligand occupations
    #   - globs: class with global variables
    # OUTPUT
    #   - strfiles: list of xyz files generated
    #   - emsg: error messages
    emsg = False
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    # get global variables class
    ############ LOAD DICTIONARIES ############
    mcores = readdict(installdir+'/Cores/cores.dict')
    licores = readdict(installdir+'/Ligands/ligands.dict')
    bindcores = readdict(installdir+'/Bind/bind.dict')
    ########## END LOAD DICTIONARIES ##########
    strfiles = []
    ########## START FUNCTIONALIZING ##########
    # load molecule core
    core,emsg = core_load(installdir,args.core,mcores)
    if emsg:
        return False,emsg
    core.convert2mol3D() # convert to mol3D
    # copy initial core for backup
    initcore3D = mol3D()
    initcore3D.copymol3D(core)
    sanity = False
    # check if ligands specified for functionalization
    if (ligands):
        # check if simple coordination complex or not
        if core.natoms == 1:
            core3D,complex3D,emsg = mcomplex(args,core,ligands,ligoc,installdir,licores,globs)
        else:
            # functionalize custom core
            core3D,emsg = customcore(args,core,ligands,ligoc,installdir,licores,globs)
        if emsg:
            return False,emsg
    else:
        core3D = initcore3D
    ############ END FUNCTIONALIZING ###########
    # generate multiple geometric arrangements
    Nogeom = int(args.bindnum) if args.bindnum and args.bind else 1 # number of different combinations
    ligname = '' # name of file
    nosmiles = 0 
    # generate name of the file
    for l in ligands:
        if l not in licores.keys():
            if '.xyz' in l or '.mol' in l:
                l = l.split('.')[-1]
                l = l.rsplit('/')[-1]
            else:
                if args.sminame:
                    if globs.nosmiles > 1:
                        ismidx = nosmiles
                    else:
                        ismidx = 0 
                    if len(args.sminame) > ismidx:
                        l = args.sminame[ismidx][0:2]
                    else:
                        l = l = 'smi'+str(nosmiles)
                else:
                    l = 'smi'+str(nosmiles)
                nosmiles += 1
        ligname += ''.join("%s" % l[0:2])
    if args.bind:
        # load bind, add hydrogens and convert to mol3D
        bind,bsmi,emsg = bind_load(installdir,args.bind,bindcores)
        if emsg:
            return False,emsg
        bind.convert2mol3D()
        an3D = bind # change name
        # get core size
        mindist = core3D.molsize()
        # assign reference point
        Rp = initcore3D.centermass()
        # Generate base case (separated structures)
        an3Db = mol3D()
        an3Db.copymol3D(an3D)
        base3D = protate(an3Db,Rp,[20*mindist,0.0,0.0])
        mols = []
        if args.bcharge:
            core3D.charge += int(args.bcharge)
        elif args.calccharge:
            core3D.charge += int(an3D.charge)
        ### check if smiles string in binding species
        if bsmi:
            if args.nambsmi: # if name specified use it in file
                fname = rootdir+'/'+core.ident[0:3]+ligname+args.nambsmi[0:2]
            else: # else use default
                fname = rootdir+'/'+core.ident[0:3]+ligname+'bsm' 
        else: # else use name from binding in dictionary
            fname = rootdir+'/'+core.ident[0:3]+ligname+bind.ident[0:2]
        # check if planar
        conats = core3D.getBondedAtomsnotH(0)
        planar,pos = False, False
        if conats > 3:
            combs = itertools.combinations(conats,4)
            for comb in combs:
                r = []
                for c in comb:
                    r.append(core3D.getAtomCoords(c))
                if checkplanar(r[0],r[1],r[2],r[3]):
                    planar = True
                    th,uax = rotation_params(r[0],r[1],r[2])
                    ueq = vecdiff(r[random.randint(0,3)],core3D.getAtomCoords(0))
                    break
        for i in range(0,Nogeom+1):        
            # generate random sequence of parameters for rotate()
            totits = 0
            while True:
                phi = random.uniform(0.0,360.0)
                theta = random.uniform(-180.0,180.0)
                if args.bphi:
                    phi = float(args.bphi)
                if args.btheta:
                    theta = float(args.btheta)
                # if specific angle is requested force angle
                if (args.place and not args.bphi and not args.btheta):
                    if ('ax' in args.place):
                        theta = 90.0
                        theta1 = -90.0
                        pos = True
                    elif ('eq' in args.place):
                        theta = 0.0
                        theta1 = 180.0
                        pos = True
                    else:
                        theta = float(args.place)
                thetax = random.uniform(0.0,360.0)
                thetay = random.uniform(0.0,360.0)
                thetaz = random.uniform(0.0,360.0)
                args.btheta = theta
                args.bphi = phi
                # translate
                an3Db = mol3D()
                an3Db.copymol3D(an3D)
                # get mask of reference atoms
                if args.bref:
                    refbP = an3D.getMask(args.bref)
                else:
                    refbP = an3D.centermass()
                if planar and pos:
                    # place axial
                    R = random.uniform(float(args.mind),float(args.maxd))
                    args.mind = R
                    args.maxd = R
                    if 'ax' in args.place:
                        newmol = setPdistanceu(an3D, refbP, core3D.getAtomCoords(0),R,uax)
                    elif 'eq' in args.place:
                        P = getPointu(core3D.getAtomCoords(0),100,ueq)
                        mindist = core3D.getfarAtomdir(P)
                        newmol = setPdistanceu(an3D, refbP, core3D.getAtomCoords(0),R+mindist,ueq)
                else:
                    # get maximum distance in the correct direction
                    Pp0 = PointTranslatetoPSph(core3D.centermass(),[0.5,0.5,0.5],[0.01,theta,phi])
                    cmcore = core3D.centermass()
                    uP = getPointu(cmcore,100,vecdiff(Pp0,cmcore)) # get far away point in space
                    mindist = core3D.getfarAtomdir(uP)
                    maxdist = mindist+float(args.maxd) # Angstrom, distance of non-interaction    
                    mindist = mindist+float(args.mind) # Angstrom, distance of non-interaction
                    R = random.uniform(mindist,maxdist) # get random distance, separated for i=0
                    # rotate and place according to distance
                    tr3D = protateref(an3Db, Rp, refbP, [R,theta,phi])
                    # rotate center of mass
                    newmol = rotateRef(tr3D,refbP,[thetax,thetay,thetaz])
                    if ('theta1' in locals()):
                        an3Db = mol3D()
                        an3Db.copymol3D(an3D)
                        tr3D2 = protateref(an3Db, Rp,refbP,[R,theta1,phi])
                        newmol2 = rotateRef(tr3D2,refbP,[thetax,thetay,thetaz])
                        d1 = tr3D.distance(core3D)
                        d2 = tr3D2.distance(core3D)
                        if (d2 > d1):
                            newmol = newmol2
                # check for overlapping
                if not(newmol.overlapcheck(core3D,1)):
                    break
                if totits > 200:
                    print "WARNING: Overlapping in molecules for file "+fname+str(i)
                    break 
                totits += 1
            if (i > 0):
                # write separate xyz file
                if args.bsep:
                    core3D.writesepxyz(newmol,fname+str(i))
                else:
                    # write new xyz file
                    newmol.writemxyz(core3D,fname+str(i))
                # append filename
                strfiles.append(fname+str(i))
                getinputargs(args,fname+str(i))
            else:
                # write new xyz file
                core3D.writexyz(fname+'R')
                # append filename
                strfiles.append(fname+'R')
                # write binding molecule file
                an3Db.writexyz(fname+'B')
                strfiles.append(fname+'B')
                del an3Db
                getinputargs(args,fname+'R')
                getinputargs(args,fname+'B')
    else:
        fname = rootdir+'/'+core.ident[0:3]+ligname
        core3D.writexyz(fname)
        strfiles.append(fname)
        getinputargs(args,fname)
    pfold = rootdir.split('/',1)[-1]
    if args.calccharge:
        args.charge = core3D.charge
    # check for molecule sanity
    sanity,d0 = core3D.sanitycheck(True)
    del core3D
    if sanity:
        print 'WARNING: Generated complex is not good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
        if args.gui:
            ssmsg = 'Generated complex in folder '+rootdir+' is no good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
            qqb = mQDialogWarn('Warning',ssmsg)
            qqb.setParent(args.gui.wmain)
    if args.gui:
        args.gui.iWtxt.setText('In folder '+pfold+' generated '+str(Nogeom)+' structures!\n'+args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    print '\nIn folder '+pfold+' generated ',Nogeom,' structures!'
    return strfiles, emsg



