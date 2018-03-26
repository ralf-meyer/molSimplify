## @file structgen.py
#  Main structure generation routine
#
#  Written by Tim Ioannidis for HJK Group
#
#  Extended by JP Janet
#
#  Revised by Terry Gani
#
#  Dpt of Chemical Engineering, MIT

from molSimplify.Scripts.geometry import *
from molSimplify.Scripts.distgeom import *
from molSimplify.Scripts.io import *
from molSimplify.Scripts.nn_prep import *
from molSimplify.Classes.globalvars import *
from molSimplify.Classes.rundiag import *
from molSimplify.Classes import globalvars
from molSimplify.Classes import mol3D
from molSimplify.Informatics.decoration_manager import*
import os, sys
from pkg_resources import resource_filename, Requirement
import openbabel, random, itertools
from numpy import log, arccos, cross, dot, pi
numpy.seterr(all='raise')

## Gets the elements in set a that are not in set b
#  @param a List with elements
#  @param b List with elements
#  @return List of elements in a that are not in b
def setdiff(a,b):
    b = set(b)
    return [aa for aa in a if aa not in b]

## Gets all possible combinations for connection atoms in geometry in the case of forced order or unknown geometry
#  @param nums List of connection atoms
#  @return List of possible backbone atom combinations
def getbackbcombsall(nums):
    bbcombs = []
    for i in range(1,len(nums)+1):
        bbcombs += list(itertools.combinations(nums,i))
    for i,tup in enumerate(bbcombs):
        bbcombs[i] = list(tup)
    return bbcombs

## Gets a combination of backbone points that satisfies denticity and updates possible combinations
#  @param backbatoms List of possible backbone atom combinations
#  @param denticity Required denticity
#  @return Selected combination, updated list of possible backbone atom combinations
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

## Gets connection atoms of SMILES string ligand
#  @param args Namespace of arguments
#  @param indsmi Index of SMILES string ligand
#  @return List of connection atoms
def getsmilescat(args,indsmi):
    tt= []  # initialize list of connection atoms
    if args.smicat and len(args.smicat)>indsmi: # get connection atom(s)
        tt = args.smicat[indsmi] # default value
    else:
        tt = [0] # default value 0 connection atom
    return tt

## Gets denticity of smiles string
#  @param args Namespace of arguments
#  @param indsmi Index of SMILES string ligand
#  @return Denticity of SMILES ligand
def getsmident(args,indsmi):
    # if denticity is specified return this
    if args.smicat and len(args.smicat) > indsmi:
        return int(len(args.smicat[indsmi]))
    # otherwise return default
    else:
        return 1

## Initializes ANN
#  @param args Namespace of arguments
#  @param ligands List of ligands
#  @param occs List of ligand occupations
#  @param dents List of ligand denticities
#  @param batslist Backbone points list
#  @param tcats List of SMILES ligand connecting atoms
#  @param licores Ligand dictionary
#  @return ANN flag, predicted BL and other attributes
def init_ANN(args,ligands,occs,dents,batslist,tcats,licores):
    # initialize ANN
    ANN_attributes = dict()
    if args.skipANN:
         print('Skipping ANN')
         ANN_flag = False
         ANN_bondl = 0
         ANN_reason = 'ANN skipped by user'
    else:
         try:
             ANN_flag,ANN_reason,ANN_attributes = ANN_preproc(args,ligands,occs,dents,batslist,tcats,licores)
             if ANN_flag:
                 ANN_bondl = ANN_attributes['ANN_bondl']
             else:
                 ANN_bondl = 0
                 if args.debug:
                     print("ANN called failed with reason: " + ANN_reason)
         except:
             print("ANN call rejected")
             ANN_reason = 'uncaught exception'
             ANN_flag = False
             ANN_bondl = 0
    return ANN_flag,ANN_bondl,ANN_reason,ANN_attributes

## Initializes core and template mol3Ds and properties
#  @param args Namespace of arguments
#  @param cpoints_required Number of connecting points required
#  @return mol3D of core, template, geometry, backbone atoms, coordination number, core reference atom index
def init_template(args,cpoints_required,globs):
    # initialize core and template
    core3D = mol3D()
    m3D = mol3D()
    # container for ordered list of core reference atoms
    corerefatoms = mol3D()
    # geometry load flag
    geom = False
    backbatoms = []
    coord = 0
    # build mode
    if args.geometry and not args.ccatoms:
        # determine geometry
        coord = int(args.coord)
        # get available geometries
        coords,geomnames,geomshorts,geomgroups = getgeoms()
        maxcoord = len(geomgroups)
        # get list of possible combinations for connecting points
        bbcombsdict = globs.bbcombs_mononuc()
        # get a default geometry
        geom = geomgroups[coord-1][0]
        # check if geometry is defined and overwrite
        if args.geometry in geomshorts:
            geom = args.geometry
        else:
            emsg = "Requested geometry not available."+"Defaulting to "+geomgroups[coord-1][0]
            if args.gui:
                qqb = mQDialogWarn('Warning',emsg)
                qqb.setParent(args.gui.wmain)
            print emsg
        # load predefined backbone coordinates
        corexyz = loadcoord(geom)
        # load backbone atom combinations
        if geom in bbcombsdict.keys() and not args.ligloc:
            backbatoms = bbcombsdict[geom]
        else:
            nums = range(1,len(corexyz))
            backbatoms = getbackbcombsall(nums)
        # distort if requested
        if args.pangles:
            corexyz = modifybackbonep(corexyz,args.pangles) # point distortion
        if args.distort:
            corexyz = distortbackbone(corexyz,args.distort) # random distortion
        # add center atom
        if args.core[0].upper()+args.core[1:] in elementsbynum:
            centeratom = args.core[0].upper()+args.core[1:]
        else:
            print('WARNING: Core is not an element. Defaulting to Fe')
            centeratom = 'Fe'
        core3D.addAtom(atom3D(centeratom,corexyz[0]))
        m3D.copymol3D(core3D)
        # add connecting points to template
        for m in range(1,coord+1):
            m3D.addAtom(atom3D('X',corexyz[m]))
            corerefatoms.addAtom(core3D.getAtom(0))
            #corerefatoms.append(0)

    # functionalize mode
    else:
        # check ccatoms
        if not args.ccatoms:
            emsg = 'Connection atoms for custom core not specified. Defaulting to 1!\n'
            print emsg
            if args.gui:
                qqb = mQDialogWarn('Warning',emsg)
                qqb.setParent(args.gui.wmain)
        ccatoms = args.ccatoms if args.ccatoms else [0]
        coord = len(ccatoms)
        if args.debug:
            print('setting ccatoms ' + str(ccatoms))

        # load core
        core,emsg = core_load(args.core)
        if emsg:
            return False,emsg
        core.convert2mol3D()
        core3D.copymol3D(core)
        m3D.copymol3D(core3D)
        for i in range(cpoints_required):
            if not args.replig:
                # not replacing ligands: add Xs to ccatoms
                # NOTE: ccatoms should be a list with # elements = cpoints_required
                cpoint = getconnection(m3D,ccatoms[i],2)
                # store core reference atom
                conatom3D = atom3D(core3D.getAtom(ccatoms[i]).sym,core3D.getAtom(ccatoms[i]).coords())
                corerefatoms.addAtom(conatom3D)
                #corerefatoms.append(ccatoms[i])
                # add connecting points to template
                m3D.addAtom(atom3D(Sym='X',xyz=cpoint))
            else:
                try:
                    # replacing ligands
                    cpoint = core3D.getAtom(ccatoms[i]).coords()
                    conatoms = core3D.getBondedAtoms(ccatoms[i])
                    # find smaller submolecule, i.e., ligand to remove
                    minmol = 10000
                    mindelats = []
                    atclose = 0
                    # loop over different connected atoms
                    for cat in conatoms:
                        # find submolecule
                        delatoms = core3D.findsubMol(ccatoms[i],cat)
                        if len(delatoms) < minmol: # check for smallest
                            mindelats = delatoms
                            minmol = len(delatoms) # size
                            atclose = cat # connection atom
                        # if same atoms in ligand get shortest distance
                        elif len(delatoms)==minmol:
                            d0 = core3D.getAtom(ccatoms[i]).distance(core3D.getAtom(cat))
                            d1 = core3D.getAtom(ccatoms[i]).distance(core3D.getAtom(mindelats[0]))
                            if d0 < d1:
                                mindelats = delatoms
                                atclose = cat
                    # store core reference atom
                    conatom3D = atom3D(core3D.getAtom(atclose).sym,core3D.getAtom(atclose).coords())
                    corerefatoms.addAtom(conatom3D)
                    #corerefatoms.append(atclose)
                    delatoms = mindelats
                    # add connecting points to template
                    m3D.addAtom(atom3D(Sym='X',xyz=cpoint))
                    # for multidentate ligands: if a submolecule contains multiple ccatoms, add all of them to the template
                    for atomidx in delatoms:
                        if atomidx in ccatoms[i+1:]:
                            # add connecting points to template
                            m3D.addAtom(atom3D(Sym='X',xyz=core3D.getAtom(atomidx).coords()))
                            ccatoms.remove(atomidx)
                            corerefatoms.addAtom(conatom3D)
                    # update remaining ccatoms according to deleted atoms
                    if len(ccatoms) > i+1:
                        for cccat in range(i+1,len(ccatoms)):
                            lshift = len([a for a in delatoms if a < ccatoms[cccat]])
                            ccatoms[cccat] -= lshift
                    # delete submolecule
                    core3D.deleteatoms(delatoms)
                    m3D.deleteatoms(delatoms)

                except IndexError:
                    pass
            nums = m3D.findAtomsbySymbol('X')
            backbatoms = getbackbcombsall(nums)
    # set charge from oxidation state if desired
    if args.calccharge:
        if args.oxstate:
            if args.oxstate in romans.keys():
                core3D.charge = int(romans[args.oxstate])
            else:
                core3D.charge = int(args.oxstate)
    return m3D,core3D,geom,backbatoms,coord,corerefatoms

## Initializes ligand 3D geometry and properties
#  @param args Namespace of arguments
#  @param lig mol3D of ligand
#  @param tcats List of SMILES ligand connecting atoms
#  @param keepHs flag for keeping H atoms on connecting atoms
#  @param i Ligand index
#  @return mol3D of ligand, flag for pi-coordination, pi-coordinating atoms
def init_ligand(args,lig,tcats,keepHs,i):
    globs = globalvars()
    rempi = False
    # if SMILES string, copy connecting atoms list to mol3D properties
    if not lig.cat and tcats[i]:
        if 'c' in tcats[i]:
            lig.cat = [lig.natoms]
        else:
            lig.cat = tcats[i]

    # change name
    lig3D = mol3D()
    lig3D.copymol3D(lig)
    # check for pi-coordinating ligand
    ligpiatoms = []
    if 'pi' in lig.cat:
        lig3Dpiatoms = mol3D()
        for k in lig.cat[:-1]:
            lig3Dpiatoms.addAtom(lig3D.getAtom(k))
            lig3Dpiatoms.addAtom(lig3D.getAtom(k))
        ligpiatoms = lig.cat[:-1]
        lig3D.addAtom(atom3D('C',lig3Dpiatoms.centermass()))
        lig.cat = [lig3D.natoms-1]
        rempi = True
    # perform FF optimization if requested (not supported for pi-coordinating ligands)
    if args.ff and 'b' in args.ffoption and not rempi:
        if 'b' in lig.ffopt.lower():
            print 'FF optimizing ligand'
            lig3D,enl = ffopt(args.ff,lig3D,lig3D.cat,0,[],False,[],100,args.debug)
    # skip hydrogen removal for pi-coordinating ligands
    if not rempi:
        # check smarts match
        if 'auto' in keepHs[i]:
            for j,catom in enumerate(lig.cat):
                match = findsmarts(lig3D.OBMol,globs.remHsmarts,catom)
                if match:
                    keepHs[i][j] = False
                else:
                    keepHs[i][j] = True
        # remove one hydrogen from each connecting atom with keepH false
        for j,cat in enumerate(lig.cat):
            Hs = lig3D.getHsbyIndex(cat)
            if len(Hs) > 0 and not keepHs[i][j]:
                if args.debug:
                    print('modifying charge down from ' + str(lig3D.charge))
                    try:
                        print('debug keepHs check, removing? ' + str(keepHs) + ' i = ' +str(i)+
                    ' , j = ' +str(j) + ' lig = ' + str(lig.coords()) + ' is keephs[i] ' + str(keepHs[i] ) +
                     ' length of keepHs list  '+ str(len(keepHs)))
                    except:
                        pass
                # check for cats indices
                if cat > Hs[0]:
                    lig.cat[j] -= 1
                lig3D.deleteatom(Hs[0])
                lig3D.charge = lig3D.charge - 1
    # Conformer search for multidentate SMILES ligands
    lig3D.convert2OBMol()
    if len(lig.cat) > 1 and tcats[i]:
        lig3D = GetConf(lig3D,lig.cat)    # check if ligand should decorated
    return lig3D,rempi,ligpiatoms

## Distorts backbone according to user specified angles
#  @param backb List with points comprising the backbone
#  @param pangles Pairs of theta/phi angles in DEGREES
#  @return List of distorted backbone points
def modifybackbonep(backb, pangles):
    for i,ll in enumerate(pangles):
        if ll:
            theta = pi*float(ll.split('/')[0])/180.0
            phi = pi*float(ll.split('/')[-1])/180.0
            backb[i+1] = PointTranslateSph(backb[0],backb[i+1],[distance(backb[0],backb[i+1]),theta,phi])
    return backb

## Randomly distorts backbone
#  @param backb List with points comprising the backbone
#  @param distort % distortion of the backbone
#  @return List of distorted backbone points
def distortbackbone(backb, distort):
    for i in range(1,len(backb)):
            theta = random.uniform(0.0,0.01*int(distort)) # *0.5
            phi = random.uniform(0.0,0.01*int(distort)*0.5) # *0.5
            backb[i] = PointTranslateSph(backb[0],backb[i],[distance(backb[0],backb[i]),theta,phi])
    return backb

##  Smart reorder ligands by denticity (-ligalign True)
#   @param args Namespace of arguments
#   @param ligs List of ligands
#   @param dentl List of ligand denticities
#   @param licores Ligand dictionary
#   return Reordered ligand indices
def smartreorderligs(args,ligs,dentl,licores):
    # reorder ligands
    globs = globalvars()
    if not args.ligalign:
        indcs = range(0,len(ligs))
        return indcs
    lsizes = []
    for ligand in ligs:
        lig,emsg = lig_load(ligand) # load ligand
        lig.convert2mol3D()
        lsizes.append(lig.natoms)
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

## Main constrained FF opt routine
#
#  To optimize metal-containing complexes with MMFF94, an intricate procedure of masking the metal atoms and manually editing their valences is applied.
#
#  OpenBabel's implementation of MMFF94 may run extremely slowly on some systems. If so, consider switching to UFF.
#
#  @param ff Force field to use, available MMFF94, UFF, Ghemical, GAFF
#  @param mol mol3D of molecule to be optimized
#  @param connected List of indices of connection atoms to metal
#  @param constopt Flag for constrained optimization - 0: unconstrained, 1: fixed connecting atom positions, 2: fixed connecting atom distances
#  @param frozenats List of frozen atom indices
#  @param frozenangles Flag for frozen angles, equivalent to constopt==1
#  @param mlbonds List of M-L bonds for distance constraints
#  @param nsteps Number of steps to take - Adaptive: run only enough steps to remove clashes, default 200
#  @param debug Flag for debug info printing
#  @return FF-calculated energy, mol3D of optimized molecule
def ffopt(ff,mol,connected,constopt,frozenats,frozenangles,mlbonds,nsteps,debug=False):
    globs = globalvars()
    metals = range(21,31)+range(39,49)+range(72,81)
    ### check requested force field
    ffav = 'mmff94, uff, ghemical, gaff, mmff94s' # force fields
    if ff.lower() not in ffav:
        print 'Requested force field not available. Defaulting to MMFF94'
        ff = 'mmff94'
    # perform constrained ff optimization if requested after #
    if (constopt > 0):
        # get metal
        midx = mol.findMetal()
        # convert mol3D to OBMol
        mol.convert2OBMol()
        OBMol = mol.OBMol
        # initialize force field
        forcefield = openbabel.OBForceField.FindForceField(ff)
        # initialize constraints
        constr = openbabel.OBFFConstraints()
        # openbabel indexing starts at 1 !!!
        # convert metals to carbons for FF
        indmtls = []
        mtlsnums = []
        for iiat,atom in enumerate(openbabel.OBMolAtomIter(OBMol)):
            if atom.GetAtomicNum() in metals:
                indmtls.append(iiat)
                mtlsnums.append(atom.GetAtomicNum())
                atom.SetAtomicNum(6)
        # freeze and ignore metals
        for midxm in indmtls:
            constr.AddAtomConstraint(midxm+1) # indexing babel
        # add coordinating atom constraints
        for ii,catom in enumerate(connected):
            if constopt==1 or frozenangles:
                constr.AddAtomConstraint(catom+1) # indexing babel
            else:
                constr.AddDistanceConstraint(midx+1,catom+1,mlbonds[ii]) # indexing babel
        bridgingatoms = []
        # identify bridging atoms in the case of bimetallic cores, as well as single-atom ligands (oxo, nitrido)
        # these are immune to deletion
        for i in range(mol.natoms):
            nbondedmetals = len([idx for idx in range(len(mol.getBondedAtoms(i))) if mol.getAtom(mol.getBondedAtoms(i)[idx]).ismetal()])
            if nbondedmetals > 1 or (nbondedmetals == 1 and len(mol.getBondedAtoms(i)) == 1):
                bridgingatoms.append(i)
        # ensure correct valences for FF setup
        for m in indmtls:
            # first delete all metal-ligand bonds excluding bridging atoms
            for i in range(len(mol.getBondedAtoms(m))):
                if OBMol.GetBond(m+1,mol.getBondedAtoms(m)[i]+1) is not None and mol.getBondedAtoms(m)[i] not in bridgingatoms:
                    OBMol.DeleteBond(OBMol.GetBond(m+1,mol.getBondedAtoms(m)[i]+1))
            # then add back one metal-ligand bond for FF
            if OBMol.GetAtom(m+1).GetValence() == 0:
                for i in mol.getBondedAtomsOct(m,1+len(bridgingatoms)):
                    if OBMol.GetAtom(m+1).GetValence() < 1 and i not in bridgingatoms:
                        OBMol.AddBond(m+1,i+1,1)
        # freeze small ligands
        for cat in frozenats:
            constr.AddAtomConstraint(cat+1) # indexing babel
        #if debug:
        #    for iiat,atom in enumerate(openbabel.OBMolAtomIter(OBMol)):
        #        print ('atom '+str(iiat)+' atomic num '+str(atom.GetAtomicNum())+' valence '+str(atom.GetValence()))
        # set up forcefield
        s = forcefield.Setup(OBMol,constr)
        if s == False:
            print('FF setup failed')
        # force field optimize structure
        elif nsteps == 'Adaptive':
            i = 0
            while i < 20:
                forcefield.ConjugateGradients(50)
                forcefield.GetCoordinates(OBMol)
                mol.OBMol = OBMol
                mol.convert2mol3D()
                overlap,mind = mol.sanitycheck(True)
                if not overlap:
                    break
                i += 1
        elif nsteps != 0:
            try:
                n = nsteps
            except:
                n = 100
            forcefield.ConjugateGradients(n)
            forcefield.GetCoordinates(OBMol)
            mol.OBMol = OBMol
            mol.convert2mol3D()
        else:
            forcefield.GetCoordinates(OBMol)
        en = forcefield.Energy()
        mol.OBMol = OBMol
        # reset atomic number to metal
        for i,iiat in enumerate(indmtls):
            mol.OBMol.GetAtomById(iiat).SetAtomicNum(mtlsnums[i])
        mol.convert2mol3D()
        del forcefield, constr, OBMol
    else:
        # initialize constraints
        constr = openbabel.OBFFConstraints()
        # add atom constraints
        for catom in connected:
            constr.AddAtomConstraint(catom+1) # indexing babel
        # set up forcefield
        forcefield = openbabel.OBForceField.FindForceField(ff)
        #if len(connected) < 2:
            #mol.OBMol.localopt('mmff94',100) # add hydrogens and coordinates
        OBMol = mol.OBMol # convert to OBMol
        s = forcefield.Setup(OBMol,constr)
        # force field optimize structure
        if OBMol.NumHvyAtoms() > 10:
            forcefield.ConjugateGradients(50)
        else:
            forcefield.ConjugateGradients(200)
        forcefield.GetCoordinates(OBMol)
        en = forcefield.Energy()
        mol.OBMol = OBMol
        mol.convert2mol3D()
        del forcefield, constr, OBMol
    return mol,en

## Finds the optimum attachment point for an atom/group to a central atom given the desired bond length
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnection(core,cidx,BL):
    ncore = core.natoms
    groups = core.getBondedAtoms(cidx)
    ccoords = core.getAtom(cidx).coords()
    # brute force search
    cpoint = []
    objopt = 0
    for itheta in range(1,359,1):
        for iphi in range(1,179,1):
            P = PointTranslateSph(ccoords,ccoords,[BL,itheta,iphi])
            dists = []
            for ig in groups:
                dists.append(distance(core.getAtomCoords(ig),P))
            obj = min(dists)
            if obj > objopt:
                objopt = obj
                cpoint = P
    return cpoint

## Checks if connecting atom of lig3D is part of SMARTS pattern
#  @param lig3D OBMol of mol3D
#  @param smarts List of SMARTS patterns (strings)
#  @param catom Connecting atom of lig3D (zero based numbering)
#  @return SMARTS match flag
def findsmarts(lig3D,smarts,catom):
    mall = []
    for smart in smarts:
        # initialize SMARTS matcher
        sm = openbabel.OBSmartsPattern()
        sm.Init(smart)
        sm.Match(lig3D)
        matches = list(sm.GetUMapList())
        # unpack tuple
        matches = [i for sub in matches for i in sub]
        for m in matches:
            if m not in mall:
                mall.append(m)
    if catom+1 in mall:
        return True
    else:
        return False

## Aligns a ligand's center of symmetry along the metal-connecting atom axis
#  @param corerefcoords Core reference coordinates
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param core3D mol3D of partially built complex
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand
def align_lig_centersym(corerefcoords,lig3D,atom0,core3D,EnableAutoLinearBend):
    # rotate to align center of symmetry
    globs = globalvars()
    r0 = corerefcoords
    r1 = lig3D.getAtom(atom0).coords()
    lig3Db = mol3D()
    lig3Db.copymol3D(lig3D)
    auxmol = mol3D()
    for at in lig3D.getBondedAtoms(atom0):
        auxmol.addAtom(lig3D.getAtom(at))
    r2 = auxmol.centersym()
    theta,u = rotation_params(r0,r1,r2)
    # rotate around axis and get both images
    lig3D = rotate_around_axis(lig3D,r1,u,theta)
    lig3Db = rotate_around_axis(lig3Db,r1,u,theta-180)
    # compare shortest distances to core reference coordinates
    d2 = distance(r0,lig3D.centersym())
    d1 = distance(r0,lig3Db.centersym())
    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
    # additional rotation for bent terminal connecting atom:
    if auxmol.natoms == 1:
        if distance(auxmol.getAtomCoords(0),lig3D.getAtomCoords(atom0)) > 0.8*(auxmol.getAtom(0).rad + lig3D.getAtom(atom0).rad) and EnableAutoLinearBend:
            print('bending of linear terminal ligand')
            ##warning: force field might overwrite this
            r1 = lig3D.getAtom(atom0).coords()
            r2 = auxmol.getAtom(0).coords()
            theta,u = rotation_params([1,1,1],r1,r2)
            lig3D = rotate_around_axis(lig3D,r1,u,globs.linearbentang)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Aligns a linear pi ligand's connecting point to the metal-ligand axis
#  @param corerefcoords Core reference coordinates
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param ligpiatoms List of ligand pi-connecting atom indices
#  @return mol3D of aligned ligand
def align_linear_pi_lig(corerefcoords,lig3D,atom0,ligpiatoms):
    # first rotate in the metal plane to ensure perpendicularity
    r0 = corerefcoords
    r1 = lig3D.getAtom(ligpiatoms[0]).coords()
    r2 = lig3D.getAtom(ligpiatoms[1]).coords()
    theta,u = rotation_params(r0,r1,r2)
    objfuncopt = 90
    thetaopt = 0
    for theta in range(0,360,1):
        lig3D_tmp = mol3D()
        lig3D_tmp.copymol3D(lig3D)
        lig3D_tmp = rotate_around_axis(lig3D_tmp, lig3D_tmp.getAtom(atom0).coords(), u, theta)
        #objfunc = abs(vecangle(vecdiff(lig3D_tmp.getAtom(atom0).coords(),corerefcoords),vecdiff(lig3D_tmp.getAtom(ligpiatoms[0]).coords(),lig3D_tmp.getAtom(ligpiatoms[1]).coords()))-90)
        objfunc = abs(distance(lig3D_tmp.getAtom(ligpiatoms[0]).coords(),corerefcoords) - distance(lig3D_tmp.getAtom(ligpiatoms[1]).coords(),corerefcoords))
        if objfunc < objfuncopt:
            thetaopt = theta
            objfuncopt = objfunc
            lig3Dopt = mol3D() # lig3Dopt = lig3D_tmp DOES NOT WORK!!!
            lig3Dopt.copymol3D(lig3D_tmp)
    lig3D = lig3Dopt
    # then rotate 90 degrees about the bond axis to further reduce steric repulsion
    r1 = lig3D.getAtom(ligpiatoms[0]).coords()
    r2 = lig3D.getAtom(ligpiatoms[1]).coords()
    u = vecdiff(r1,r2)
    lig3D_tmpa = mol3D()
    lig3D_tmpa.copymol3D(lig3D)
    lig3D_tmpa = rotate_around_axis(lig3D_tmpa, lig3D_tmpa.getAtom(atom0).coords(), u, 90)
    lig3D_tmpb = mol3D()
    lig3D_tmpb.copymol3D(lig3D)
    lig3D_tmpb = rotate_around_axis(lig3D_tmpb, lig3D_tmpb.getAtom(atom0).coords(), u, -90)
    d1 = distance(corerefcoords,lig3D_tmpa.centermass())
    d2 = distance(corerefcoords,lig3D_tmpb.centermass())
    #lig3D = lig3D if (d1 < d2)  else lig3Db
    lig3D = lig3D_tmpa if (d1 > d2) else lig3D_tmpb # pick the better structure
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Checks if ligand has a linear coordination environment (e.g., OCO) and ensures perpendicularity to M-L axis
#  @param corerefcoords Core reference coordinates
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @return mol3D of rotated ligand
def check_rotate_linear_lig(corerefcoords,lig3D,atom0):
    auxm = mol3D()
    lig3D_aligned = mol3D()
    for at in lig3D.getBondedAtoms(atom0):
        auxm.addAtom(lig3D.getAtom(at))
    if auxm.natoms > 1:
        r0 = lig3D.getAtom(atom0).coords()
        r1 = auxm.getAtom(0).coords()
        r2 = auxm.getAtom(1).coords()
        if checkcolinear(r1,r0,r2):
            # rotate so that O-C-O bond is perpendicular to M-L axis
            theta,urot = rotation_params(r1,corerefcoords,r2)
            theta = vecangle(vecdiff(r0,corerefcoords),urot)
            lig3D = rotate_around_axis(lig3D,r0,urot,theta)
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Checks if ligand has is symmetric about connecting atom (center of symmetry coincides with connecting atom) and minimizes clashes with rest of complex
#  @param corerefcoords Core reference coordinates
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param core3D mol3D of partially built complex
#  @return mol3D of rotated ligand
def check_rotate_symm_lig(corerefcoords,lig3D,atom0,core3D):
    if distance(lig3D.getAtom(atom0).coords(),lig3D.centersym()) < 8.0e-2:
        at = lig3D.getBondedAtoms(atom0)
        r0 = lig3D.getAtom(atom0).coords()
        r1 = lig3D.getAtom(at[0]).coords()
        r2 = lig3D.getAtom(at[1]).coords()
        theta,u = rotation_params(r0,r1,r2)
        theta = vecangle(u,vecdiff(r0,corerefcoords))
        urot = cross(u,vecdiff(r0,corerefcoords))
        # rotate around axis and get both images
        lig3Db = mol3D()
        lig3Db.copymol3D(lig3D)
        lig3D = rotate_around_axis(lig3D,r0,urot,theta)
        lig3Db = rotate_around_axis(lig3Db,r0,urot,-theta)
        # compute shortest distances to core
        d2 = lig3D.mindist(core3D)
        d1 = lig3Db.mindist(core3D)
        lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Rotates aligned ligand about M-L axis to minimize steric clashes with rest of complex
#  @param corerefcoords Core reference coordinates
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param core3D mol3D of partially built complex
#  @return mol3D of rotated ligand
def rotate_MLaxis_minimize_steric(corerefcoords,lig3D,atom0,core3D):
    r1 = lig3D.getAtom(atom0).coords()
    u = vecdiff(r1,corerefcoords)
    dtheta = 2
    optmax = -9999
    totiters = 0
    lig3Db = mol3D()
    lig3Db.copymol3D(lig3D)
    # maximize a combination of minimum distance between atoms and center of mass distance
    while totiters < 180:
        lig3D = rotate_around_axis(lig3D,r1,u,dtheta)
        d0 = lig3D.mindist(core3D) # shortest distance
        d0cm = lig3D.distance(core3D) # center of mass distance
        iteropt = d0cm+10*log(d0)
        if (iteropt > optmax): # if better conformation, keep
            lig3Db = mol3D()
            lig3Db.copymol3D(lig3D)
            optmax = iteropt
        totiters += 1
    lig3D = lig3Db
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Rotates a connecting atom of a multidentate ligand to improve H atom placement
#
#  There are separate routines for terminal connecting atoms and intermediate connecting atoms.
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @param n Index of connecting atom
#  @param mcoords Core reference (usually a metal) coordintes
#  @param core3D mol3D of partially built complex
#  @return mol3D of rotated ligand
def rotate_catom_fix_Hs(lig3D,catoms,n,mcoords,core3D):
    # isolate fragment to be rotated
    confrag3D = mol3D()
    confragatomlist = []
    danglinggroup = []
    catoms_other = catoms[:]
    catoms_other.pop(n)
    # add connecting atom
    confrag3D.addAtom(lig3D.getAtom(catoms[n]))
    confragatomlist.append(catoms[n])
    # add all Hs bound to connecting atom
    for ii in lig3D.getHsbyIndex(catoms[n]):
        confrag3D.addAtom(lig3D.getAtom(ii))
        confragatomlist.append(ii)
    # add dangling groups
    anchoratoms = []
    for atom in lig3D.getBondedAtomsnotH(catoms[n]):
        subm = lig3D.findsubMol(atom,catoms[n])
        if len(list(set(subm).intersection(catoms_other))) == 0:
            danglinggroup = subm
        else:
            bridginggroup = subm
            if list(set(subm).intersection(lig3D.getBondedAtoms(catoms[n])))[0] not in anchoratoms:
                anchoratoms.append(list(set(subm).intersection(lig3D.getBondedAtoms(catoms[n])))[0])
    for atom in danglinggroup:
        confrag3D.addAtom(lig3D.getAtom(atom))
        confragatomlist.append(atom)
    if confrag3D.natoms > 1:
        # terminal connecting atom
        confrag3Dtmp = mol3D()
        confrag3Dtmp.copymol3D(confrag3D)
        if len(anchoratoms) == 1:
            anchoratom = anchoratoms[0]
            anchor = lig3D.getAtomCoords(anchoratom)
            if not checkcolinear(anchor,confrag3D.getAtomCoords(0),confrag3D.getAtomCoords(1)):
                refpt = confrag3D.getAtomCoords(0)
                u = vecdiff(refpt,anchor)
                dtheta = 5
                objs = []
                objopt = 0
                localmaxs = []
                thetas = range(0,360,dtheta)
                for theta in thetas:
                    confrag3Dtmp = rotate_around_axis(confrag3Dtmp,refpt,u,dtheta)
                    auxmol1 = mol3D()
                    auxmol1.addAtom(confrag3Dtmp.getAtom(0))
                    for at in confrag3Dtmp.getBondedAtoms(0):
                        auxmol1.addAtom(confrag3Dtmp.getAtom(at))
                    auxmol1.addAtom(lig3D.getAtom(anchoratom))
                    auxmol2 = mol3D()
                    auxmol2.copymol3D(confrag3Dtmp)
                    #objs.append(distance(mcoords,auxmol.centersym()))
                    if auxmol2.natoms > 3:
                        obj = auxmol2.mindisttopoint(mcoords)
                    else:
                        obj = distance(mcoords,auxmol1.centersym())
                    if obj > objopt:
                        objopt = obj
                        thetaopt = theta
                #for i,obj in enumerate(objs):
                    #try:
                        #if objs[i] > objs[i-1] and objs[i] > objs[i+1]:
                            #localmaxs.append(thetas[i])
                    #except IndexError:
                        #pass
            ## in future, compare multiple local maxima
            #if localmaxs == []:
                #localmaxs = [0]
            confrag3D = rotate_around_axis(confrag3D,refpt,u,thetaopt)
            #confrag3D = rotate_around_axis(confrag3D,refpt,u,localmaxs[0])
        # non-terminal connecting atom
        elif len(anchoratoms) == 2:
            refpt = confrag3D.getAtomCoords(0)
            anchorcoords1 = lig3D.getAtomCoords(anchoratoms[0])
            anchorcoords2 = lig3D.getAtomCoords(anchoratoms[1])
            u = vecdiff(anchorcoords1,anchorcoords2)
            dtheta = 5
            objs = []
            localmaxs = []
            thetas = range(0,360,dtheta)
            for theta in thetas:
                confrag3Dtmp = rotate_around_axis(confrag3Dtmp,refpt,u,dtheta)
                newHcoords = confrag3Dtmp.getAtomCoords(1)
                objs.append(distance(newHcoords,anchorcoords1)+distance(newHcoords,anchorcoords2)+distance(newHcoords,mcoords))
            for i,obj in enumerate(objs):
                try:
                    if objs[i] > objs[i-1] and objs[i] > objs[i+1]:
                        localmaxs.append(thetas[i])
                except IndexError:
                    pass
            if localmaxs == []:
                localmaxs = [0]
            confrag3D = rotate_around_axis(confrag3D,refpt,u,localmaxs[0])
        for i,atom in enumerate(confragatomlist):
            lig3D.getAtom(atom).setcoords(confrag3D.getAtomCoords(i))
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Rotates connecting atoms of multidentate ligands to improve H atom placement
#
#  Loops over rotate_catom_fix_Hs().
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @param mcoords Core reference (usually a metal) coordintes
#  @param core3D mol3D of partially built complex
#  @return mol3D of rotated ligand
def rotate_catoms_fix_Hs(lig3D,catoms,mcoords,core3D):
    for i,n in enumerate(catoms):
        #if len(lig3D.getHsbyIndex(n)) > 0:
        lig3D = rotate_catom_fix_Hs(lig3D,catoms,i,mcoords,core3D)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned

## Finds and rotates a rotatable bond in a bidentate ligand to obtain cis conformer
#
#  For certain ligands (e.g., NCCN) this is a more efficient version of the generic distance geometry conformer search.
#
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @return mol3D of rotated ligand
def find_rotate_rotatable_bond(lig3D,catoms):
    bats = list(set(lig3D.getBondedAtomsnotH(catoms[0])) | set(lig3D.getBondedAtomsnotH(catoms[1])))
    rb1 = 1000
    rb2 = 1000
    for ii in range(lig3D.OBMol.NumBonds()):
        bd = lig3D.OBMol.GetBond(ii)
        bst = bd.GetBeginAtomIdx()
        ben = bd.GetEndAtomIdx()
        if bd.IsRotor() and (bst-1 in bats) and (ben-1 in bats):
            print('Rotatable bond found')
            rb1 = bst-1
            rb2 = ben-1
            break
    if (rb1 != 1000) and (rb2 != 1000): # rotatable bond present, execute rotations
        rotfrag3D = mol3D()
        # create submolecule containing atoms to be rotated (the one containing catoms[0] which is aligned first)
        subm1 = lig3D.findsubMol(rb1,rb2)
        subm2 = lig3D.findsubMol(rb2,rb1)
        if catoms[1] in subm1:
            subm = subm1
            anchor = lig3D.getAtomCoords(rb2)
            refpt = lig3D.getAtomCoords(rb1)
        elif catoms[0] in subm1:
            subm = subm2
            anchor = lig3D.getAtomCoords(rb1)
            refpt = lig3D.getAtomCoords(rb2)
        ncoord = 0
        for nii,ii in enumerate(subm):
            rotfrag3D.addAtom(lig3D.getAtom(ii))
            # find coordinating atom in submolecule
            if ii in catoms:
                ncoord = nii
        u = vecdiff(refpt,anchor)
        dtheta = 10
        theta = 0
        thetaopt = 0
        objopt = 1000
        while theta < 360: # minimize distance between connecting atoms
            rotfrag3D = rotate_around_axis(rotfrag3D,anchor,u,dtheta)
            obj = distance(lig3D.getAtomCoords(catoms[1]),rotfrag3D.getAtomCoords(ncoord))
            obj = obj + distance(lig3D.getAtomCoords(catoms[0]),rotfrag3D.getAtomCoords(ncoord))
            if obj < objopt:
                thetaopt = theta
                objopt = obj
            theta = theta + dtheta
        rotfrag3D = rotate_around_axis(rotfrag3D,anchor,u,thetaopt)
        jj = 0
        for ii in subm: # replace coordinates
            lig3D.getAtom(ii).setcoords(rotfrag3D.getAtomCoords(jj))
            jj = jj + 1
    return lig3D

## Gets target M-L distance from desired source (custom, sum cov rad or ANN)
## Aligns a monodentate ligand to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param ligand Name of ligand for dictionary lookup
#  @param metal atom3D of atom 1 (usually a metal)
#  @param MLb Custom M-L bond length (if any)
#  @param i Ligand serial number
#  @param ANN_flag Flag for ANN activation
#  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @return M-L bond length in Angstroms
def get_MLdist(args,lig3D,atom0,ligand,metal,MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds):
    # first check for user-specified distances and use them
    if MLb and MLb[i]:
        print('using user-specified M-L distances')
        if 'c' in MLb[i].lower():
            bondl = metal.rad + lig3D.getAtom(atom0).rad
        else:
            bondl = float(MLb[i])
    else:
    # otherwise, check for exact DB match
        bondl,exact_match = get_MLdist_database(args,metal,lig3D,atom0,ligand,MLbonds)
        try:
            this_diag.set_dict_bl(bondl)
        except:
            pass
        if not exact_match and ANN_flag:
            # if no exact match found and ANN enabled, use it
            print('no M-L match in DB, using ANN')
            bondl =  ANN_bondl
        elif exact_match:
            print('using exact M-L match from DB')
        else:
            print('Warning: ANN not active and exact M-L match not found in DB, distance may not be accurate')
    return bondl

## Loads M-L bond length from database and reports if compound is in DB
#  @param args Namespace of arguments
#  @param metal atom3D of atom 1 (usually a metal)
#  @param lig3D mol3D of ligand
#  @param atom0 Ligand connecting atom index
#  @param ligand Name of ligand
#  @param MLbonds M-L dictionary
#  @return Bond length in Angstroms, flag for exact DB match
def get_MLdist_database(args,metal,lig3D,atom0,ligand,MLbonds):
    # check for roman letters in oxstate
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
    key.append((metal.sym,oxs,spin,lig3D.getAtom(atom0).sym,ligand))
    key.append((metal.sym,oxs,spin,lig3D.getAtom(atom0).sym,'-')) # disregard exact ligand
    key.append((metal.sym,'-','-',lig3D.getAtom(atom0).sym,ligand)) # disregard oxstate/spin
    key.append((metal.sym,'-','-',lig3D.getAtom(atom0).sym,'-')) # else just consider bonding atom
    found = False
    exact_match = False
    # search for data
    for kk in key:
        if (kk in MLbonds.keys()): # if exact key in dictionary
            bondl = float(MLbonds[kk])
            found = True
            if (kk == ((metal.sym,oxs,spin,lig3D.getAtom(atom0).sym,ligand))): ## exact match
               exact_match = True
            break
    if not found: # last resort covalent radii
        bondl = metal.rad + lig3D.getAtom(atom0).rad
    if args.debug:
        print('ms default distance is  ' + str(bondl))
    return bondl,exact_match

## Get backbone atoms from template
#  @param args Namespace of arguments
#  @param batslist List of backbone connecting atoms for each ligand
#  @param ligsused Number of ligands placed
#  @return Backbone connecting atoms for ligand
def get_batoms(args,batslist,ligsused):
    batoms = batslist[ligsused]
    if len(batoms) < 1 :
        emsg = 'Connecting all ligands is not possible. Check your input!'
        if args.gui:
            qqb = mQDialogWarn('Warning',emsg)
            qqb.setParent(args.gui.wmain)
    return batoms

## Crude rotations to improve alignment of the 2nd connecting atom of a bidentate ligand
#  @param args Namespace of arguments
#  @param lig3D mol3D of ligand
#  @param core3D mol3D of partially built complex
#  @param catoms List of ligand connecting atom indices
#  @param r1 Coordinates of ligand first connecting atom
#  @param r0 Coordinates of core reference point
#  @param m3D mol3D of backbone template
#  @param batoms List of backbone atom indices
#  @param corerefcoords Coordinates of core reference atom
#  @return mol3D of aligned ligand, coordinates of second backbone point
def align_dent2_catom2_coarse(args,lig3D,core3D,catoms,r1,r0,m3D,batoms,corerefcoords):
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
    md = min(distance(r0l,corerefcoords),distance(r1l,corerefcoords))
    if lig3D.mindist(core3D) < md:
        lig3D = rotate_around_axis(lig3D,r0l,vecdiff(r1l,r0l),180.0)
    # correct plane
    r0b = m3D.getAtom(batoms[0]).coords()
    r1b = m3D.getAtom(batoms[1]).coords()
    r0l = lig3D.getAtom(catoms[0]).coords()
    r1l = lig3D.getAtom(catoms[1]).coords()
    rm = lig3D.centermass()
    urot = vecdiff(r1l,r0l)
    #theta,ub = rotation_params(corerefcoords,r0b,r1b)
    #theta,ul = rotation_params(rm,r0l,r1l)
    #if (norm(ub)*norm(ul)) > 1e-8:
        #theta = 180*arccos(dot(ub,ul)/(norm(ub)*norm(ul)))/pi-180.0
    #else:
        #theta = 0.0
    # rotate around axis
    objopt = 0
    for theta in range(0,360,5):
        lig3D_tmp = mol3D()
        lig3D_tmp.copymol3D(lig3D)
        lig3D_tmp = rotate_around_axis(lig3D_tmp,r1,urot,theta)
        lig3D_tmp2 = mol3D()
        lig3D_tmp2.copymol3D(lig3D_tmp)
        H1 = lig3D_tmp2.getBondedAtomsH(catoms[1])
        H2 = lig3D_tmp2.getBondedAtomsH(catoms[0])
        lig3D_tmp2.deleteatoms([catoms[1]]+[catoms[0]]+H1+H2)
        obj = lig3D_tmp2.mindisttopoint(corerefcoords)
        if obj > objopt:
            objopt = obj
            lig3Dopt = mol3D()
            lig3Dopt.copymol3D(lig3D_tmp)
    lig3D = mol3D()
    lig3D.copymol3D(lig3Dopt)
    tmp3D = mol3D()
    tmp3D.copymol3D(m3D)
    tmp3D.combine(lig3D)
    tmp3D.writexyz('new')
    #lig3Db = mol3D()
    #lig3Db.copymol3D(lig3D)
    #lig3D = rotate_around_axis(lig3D,r1,urot,theta)
    #lig3Db = rotate_around_axis(lig3Db,r1,urot,-theta)
    # select best
    try:
        rm0,rm1 = lig3D.centermass(),lig3Db.centermass()
        theta,ul0 = rotation_params(rm0,r0l,r1l)
        theta,ul1 = rotation_params(rm1,r0l,r1l)
        th0 = 180*arccos(dot(ub,ul0)/(norm(ub)*norm(ul0)))/pi
        th0 = min(abs(th0),abs(180-th0))
        th1 = 180*arccos(dot(ub,ul1)/(norm(ub)*norm(ul1)))/pi
        th1 = min(abs(th1),abs(180-th1))
        lig3D = lig3D if th0 < th1 else lig3Db
    except:
        pass
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned,r1b

## Aligns second connecting atom of a bidentate ligand to balance ligand strain and the desired coordination environment.
#  @param args Namespace of arguments
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @param bondl Target M-L bond length
#  @param r1 Coordinates of ligand first connecting atom
#  @param r0 Coordinates of core reference point
#  @param core3D mol3D of partially built complex
#  @param rtarget Coordinates of target point for second connecting atom
#  @param coreref atom3D of core reference atom
#  @param MLoptbds List of final M-L bond lengths
#  @return mol3D of aligned ligand
def align_dent2_catom2_refined(args,lig3D,catoms,bondl,r1,r0,core3D,rtarget,coreref,MLoptbds):
    # compute starting ligand FF energy for later comparison
    corerefcoords = coreref.coords()
    dr = vecdiff(rtarget,lig3D.getAtom(catoms[1]).coords())
    cutoff = 5 # energy threshold for ligand strain, kcal/mol
    lig3Dtmp = mol3D()
    lig3Dtmp.copymol3D(lig3D)
    lig3Dtmp,en_start = ffopt(args.ff,lig3Dtmp,[],1,[],False,[],200,args.debug)
    # take steps between current ligand position and ideal position on backbone
    nsteps = 20
    ddr = [di/nsteps for di in dr]
    ens=[]
    finished = False
    relax = False
    while True:
        lig3Dtmp = mol3D()
        lig3Dtmp.copymol3D(lig3D)
        for ii in range(0,nsteps):
            lig3Dtmp,enl = ffopt(args.ff,lig3Dtmp,[],1,[catoms[0],catoms[1]],False,[],'Adaptive',args.debug)
            ens.append(enl)
            lig3Dtmp.getAtom(catoms[1]).translate(ddr)
            # once the ligand strain energy becomes too high, stop and accept ligand position
            # or if the ideal coordinating point is reached without reaching the strain energy cutoff, stop
            if (ens[-1] - ens[0] > cutoff) or (ii == nsteps-1):
                r0,r1 = lig3Dtmp.getAtomCoords(catoms[0]),lig3Dtmp.getAtomCoords(catoms[1])
                r01 = distance(r0,r1)
                try:
                # but if ligand still cannot be aligned, instead force alignment with a huge cutoff and then relax later
                    theta1 = 180*arccos(0.5*r01/bondl)/pi
                except:
                    print('Forcing alignment...')
                    cutoff += 5000000
                    relax = True
                    break
                theta2 = vecangle(vecdiff(r1,r0),vecdiff(corerefcoords,r0))
                dtheta = theta2-theta1
                theta,urot = rotation_params(corerefcoords,r0,r1)
                lig3Dtmp = rotate_around_axis(lig3Dtmp,r0,urot,-dtheta) # rotate so that it matches bond
                finished = True
                break
        if finished:
            break
    # for long linear ligand chains, this procedure might produce the wrong ligand curvature. If so, reflect about M-L plane
    lig3Dtmpb = mol3D()
    lig3Dtmpb.copymol3D(lig3Dtmp)
    lig3Dtmpb = reflect_through_plane(lig3Dtmpb,vecdiff(midpt(lig3Dtmpb.getAtom(catoms[0]).coords(),lig3Dtmpb.getAtom(catoms[1]).coords()),corerefcoords),lig3Dtmpb.getAtom(catoms[0]).coords())
    lig3Dtmp = lig3Dtmpb if lig3Dtmp.mindist(core3D) < lig3Dtmpb.mindist(core3D) else lig3Dtmp
    if relax:
        # Relax the ligand
        lig3Dtmp,enl = ffopt(args.ff,lig3Dtmp,[catoms[1]],2,[catoms[0]],False,MLoptbds[-2:-1],200,args.debug)
        lig3Dtmp.deleteatom(lig3Dtmp.natoms-1)
    lig3Dtmp,en_final = ffopt(args.ff,lig3Dtmp,[],1,[],False,[],0,args.debug)
    if en_final - en_start > 20:
        print 'Warning: Complex may be strained. Ligand strain energy (kcal/mol) = ' + str(en_final - en_start)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3Dtmp)
    return lig3D_aligned

## Aligns a monodentate ligand to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param ligand Name of ligand for dictionary lookup
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @param rempi Flag for pi-coordinating ligand
#  @param ligpiatoms List of pi-coordinating atom indices in ligand
#  @param MLb Custom M-L bond length (if any)
#  @param ANN_flag Flag for ANN activation
#  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param i Ligand serial number
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand, updated list of M-L bond lengths
def align_dent1_lig(args,cpoint,core3D,coreref,ligand,lig3D,catoms,rempi=False,ligpiatoms=[],MLb=[],ANN_flag=False,ANN_bondl=[],this_diag=0,MLbonds=dict(),MLoptbds=[],i=0,EnableAutoLinearBend=True):
    corerefcoords = coreref.coords()
    # connection atom in lig3D
    atom0 = catoms[0]
    # translate ligand to overlap with backbone connecting point
    lig3D.alignmol(lig3D.getAtom(atom0),cpoint)
    # determine bond length (database/cov rad/ANN)
    bondl = get_MLdist(args,lig3D,atom0,ligand,coreref,MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
    MLoptbds.append(bondl)
    # align ligand to correct M-L distance
    u = vecdiff(cpoint.coords(),corerefcoords)
    lig3D = aligntoaxis2(lig3D, cpoint.coords(), corerefcoords, u, bondl)
    if rempi and len(ligpiatoms) == 2:
        # align linear (non-arom.) pi-coordinating ligand
        lig3D = align_linear_pi_lig(corerefcoords,lig3D,atom0,ligpiatoms)
    elif lig3D.natoms > 1:
        # align ligand center of symmetry
        lig3D = align_lig_centersym(corerefcoords,lig3D,atom0,core3D,EnableAutoLinearBend)
        if lig3D.natoms > 2:
            # check for linear molecule and align
            lig3D = check_rotate_linear_lig(corerefcoords,lig3D,atom0)
            # check for symmetric molecule
            lig3D = check_rotate_symm_lig(corerefcoords,lig3D,atom0,core3D)
        # rotate around M-L axis to minimize steric repulsion
        lig3D = rotate_MLaxis_minimize_steric(corerefcoords,lig3D,atom0,core3D)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned,MLoptbds

## Aligns a bidentate ligand to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param batoms List of backbone atom indices
#  @param m3D mol3D of backbone template
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param ligand Name of ligand for dictionary lookup
#  @param lig3D mol3D of ligand
#  @param catoms List of ligand connecting atom indices
#  @param MLb Custom M-L bond length (if any)
#  @param ANN_flag Flag for ANN activation
#  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param frozenats Atoms frozen in FF optimization
#  @param i Ligand serial number
#  @return mol3D of aligned ligand, updated lists of frozen atoms and M-L bond lengths
def align_dent2_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl,this_diag,MLbonds,MLoptbds,frozenats,i):
    corerefcoords = coreref.coords()
    r0 = corerefcoords
    # get cis conformer by rotating rotatable bonds
    #lig3D = find_rotate_rotatable_bond(lig3D,catoms)
    # connection atom
    atom0 = catoms[0]
    # translate ligand to match first connecting atom to backbone connecting point
    lig3D.alignmol(lig3D.getAtom(atom0),cpoint)
    r1 = lig3D.getAtom(atom0).coords()
    # Crude rotations to bring the 2nd connecting atom closer to its ideal location
    lig3D,r1b = align_dent2_catom2_coarse(args,lig3D,core3D,catoms,r1,r0,m3D,batoms,corerefcoords)
    ## get bond length
    bondl = get_MLdist(args,lig3D,atom0,ligand,coreref,MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
    MLoptbds.append(bondl)
    MLoptbds.append(bondl)
    lig3D,dxyz = setPdistance(lig3D, r1, r0, bondl)
    # get target point for 2nd connecting atom
    rtarget = getPointu(corerefcoords, bondl, vecdiff(r1b,corerefcoords)) # get second point target
    if args.ff:
        # align 2nd connecting atom while balancing the desired location and ligand strain
        lig3D = align_dent2_catom2_refined(args,lig3D,catoms,bondl,r1,r0,core3D,rtarget,coreref,MLoptbds)
    else:
        print 'Warning: Ligand FF optimization is inactive.'
    # rotate connecting atoms to align Hs properly
    lig3D = rotate_catoms_fix_Hs(lig3D,catoms,corerefcoords,core3D)
    # freeze local geometry
    lats = lig3D.getBondedAtoms(catoms[0])+lig3D.getBondedAtoms(catoms[1])
    for lat in list(set(lats)):
        frozenats.append(lat+core3D.natoms)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned,frozenats,MLoptbds

def align_dent3_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl,this_diag,MLbonds,MLoptbds,frozenats,i):
    atom0 = catoms[1]
    corerefcoords = coreref.coords()
    # align molecule according to connection atom and shadow atom
    lig3D.alignmol(lig3D.getAtom(atom0),m3D.getAtom(batoms[1]))
    # 1. align ligand connection atoms center of symmetry
    auxm = mol3D()
    auxm.addAtom(lig3D.getAtom(catoms[0]))
    auxm.addAtom(lig3D.getAtom(catoms[2]))
    r0 = core3D.getAtom(0).coords()
    lig3Db = mol3D()
    lig3Db.copymol3D(lig3D)
    theta,urot = rotation_params(r0,lig3D.getAtom(atom0).coords(),auxm.centersym())
    lig3D = rotate_around_axis(lig3D,lig3D.getAtom(atom0).coords(),urot,theta)
    # 2. align with correct plane
    rl0,rl1,rl2 = lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),lig3D.getAtom(catoms[2]).coords()
    rc0,rc1,rc2 = m3D.getAtom(batoms[0]).coords(),m3D.getAtom(batoms[1]).coords(),m3D.getAtom(batoms[2]).coords()
    theta0,ul = rotation_params(rl0,rl1,rl2)
    theta1,uc = rotation_params(rc0,rc1,rc2)
    urot = vecdiff(rl1,corerefcoords)
    theta = vecangle(ul,uc)
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
    # 3. correct if not symmetric
    theta0,urotaux = rotation_params(lig3D.getAtom(catoms[0]).coords(),lig3D.getAtom(catoms[1]).coords(),core3D.getAtom(0).coords())
    theta1,urotaux = rotation_params(lig3D.getAtom(catoms[2]).coords(),lig3D.getAtom(catoms[1]).coords(),core3D.getAtom(0).coords())
    dtheta = 0.5*(theta1-theta0)
    if abs(dtheta) > 0.5:
        lig3D = rotate_around_axis(lig3D,lig3D.getAtom(atom0).coords(),urot,dtheta)
    # 4. flip for correct stereochemistry
    urot = vecdiff(lig3D.getAtom(catoms[1]).coords(),core3D.getAtom(0).coords())
    lig3Db = mol3D()
    lig3Db.copymol3D(lig3D)
    lig3Db = rotate_around_axis(lig3Db,lig3Db.getAtom(catoms[1]).coords(),urot,180)
    d1 = min(distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords()),distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[0]).coords()))
    d2 = min(distance(lig3Db.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords()),distance(lig3Db.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[0]).coords()))
    lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
    # 5. flip to align 1st and 3rd connection atoms
    lig3Db = mol3D()
    lig3Db.copymol3D(lig3D)
    theta,urot = rotation_params(lig3Db.getAtom(catoms[0]).coords(),lig3Db.getAtom(catoms[1]).coords(),lig3Db.getAtom(catoms[2]).coords())
    lig3Db = rotate_around_axis(lig3Db,lig3Db.getAtom(catoms[1]).coords(),urot,180)
    d1 = min(distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords()),distance(lig3D.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[0]).coords()))
    d2 = min(distance(lig3Db.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[2]).coords()),distance(lig3Db.getAtom(catoms[2]).coords(),m3D.getAtom(batoms[0]).coords()))
    lig3D = lig3D if d1 < d2 else lig3Db
    bondl = get_MLdist(args,lig3D,atom0,ligand,m3D.getAtom(0),MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
    for iib in range(0,3):
        MLoptbds.append(bondl)
    # set correct distance
    setPdistance(lig3D, lig3D.getAtom(atom0).coords(), m3D.getAtom(0).coords(), bondl)
    # rotate connecting atoms to align Hs properly
    lig3D = rotate_catoms_fix_Hs(lig3D,[catoms[0],catoms[1],catoms[2]],m3D.getAtom(0).coords(),core3D)
    # freeze local geometry
    lats = lig3D.getBondedAtoms(catoms[0])+lig3D.getBondedAtoms(catoms[1])
    for lat in list(set(lats)):
        frozenats.append(lat+core3D.natoms)
    lig3D_aligned = mol3D()
    lig3D_aligned.copymol3D(lig3D)
    return lig3D_aligned,frozenats,MLoptbds    
    
## Test function, not used anywhere
#  Tries to break a tridentate ligand at the middle connecting atom into the middle fragment and two bidentate ligands which are then aligned separately.
#
#  Will fail if there are no suitable rotatable bonds, such as in a polycyclic aromatic ligand, and return a flag indicating as such.
#  @param lig3D mol3D of tridentate ligand
#  @param catoms List of connecting atoms
#  @return mol3Ds of fragments, connecting atoms, success flag
def cleave_tridentate(lig3D,catoms):
    # end fragments
    lig3D1 = mol3D()
    catoms1 = []
    lig3D2 = mol3D()
    catoms2 = []
    # middle fragment
    lig3D3 = mol3D()
    lig3D3.addAtom(lig3D.getAtom(catoms[1]))
    status = False
    for i in lig3D.getBondedAtomsnotH(catoms[1]):
        if catoms[0] in lig3D.findsubMol(i,catoms[1]):
            catoms1.append(lig3D.findsubMol(i,catoms[1]).index(catoms[0]))
            for n,j in enumerate(lig3D.findsubMol(i,catoms[1])):
                lig3D1.addAtom(lig3D.getAtom(j))
            lig3D1.addAtom(lig3D.getAtom(catoms[1]))
            catoms1.append(lig3D1.natoms-1)
            status = True
        elif catoms[2] in lig3D.findsubMol(i,catoms[1]):
            catoms2.append(lig3D.findsubMol(i,catoms[1]).index(catoms[2]))
            for j in lig3D.findsubMol(i,catoms[1]):
                lig3D2.addAtom(lig3D.getAtom(j))
            lig3D2.addAtom(lig3D.getAtom(catoms[1]))
            catoms2.append(lig3D2.natoms-1)
            status = True
        else:
            lig3D3.addAtom(lig3D.getAtom(i))
    if lig3D1.natoms >= lig3D.natoms or lig3D2.natoms >= lig3D.natoms:
        status = False
    catoms1.reverse()
    catoms2.reverse()
    return lig3D1,catoms1,lig3D2,catoms2,lig3D3,status

## Main ligand placement routine
#  @param args Namespace of arguments
#  @param ligs List of ligands
#  @param ligoc List of ligand occupations
#  @param licores Ligand dictionary
#  @param globs Global variables
#  @return mol3D of built complex, list of all mol3D ligands and core, error messages
def mcomplex(args,ligs,ligoc,licores,globs):
    this_diag = run_diag()
    if globs.debug:
        print '\nGenerating complex with ligands and occupations:',ligs,ligoc
    if args.gui:
        args.gui.iWtxt.setText('\nGenerating complex with core:'+args.core+' and ligands: '+ ' '.join(ligs)+'\n'+args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    # initialize variables
    emsg, complex3D = False, []
    occs0 = []      # occurrences of each ligand
    toccs = 0       # total occurrence count (number of ligands)
    catsmi = []     # SMILES ligands connection atoms
    smilesligs = 0  # count how many smiles strings
    cats0 = []      # connection atoms for ligands
    dentl = []      # denticity of ligands
    connected = []  # indices in core3D of ligand atoms connected to metal
    frozenats = []  # atoms to be frozen in optimization
    freezeangles = False # custom angles imposed
    MLoptbds = []   # list of bond lengths
    rempi = False   # remove dummy pi orbital center of mass atom
    backbatoms = []
    batslist = []
    bats = []
    # load bond data
    MLbonds = loaddata('/Data/ML.dat')
    # calculate occurrences, denticities etc for all ligands
    for i,ligname in enumerate(ligs):
        # if not in cores -> smiles/file
        if ligname not in licores.keys():
            if args.smicat and len(args.smicat)>= (smilesligs+1):
                if 'pi' in args.smicat[smilesligs]:
                    cats0.append(['c'])
                else:
                    cats0.append(args.smicat[smilesligs])
            else:
                cats0.append([0])
            dent_i = len(cats0[-1])
            smilesligs += 1
        else:
            cats0.append(False)
        # otherwise get denticity from ligands dictionary
            if 'pi' in licores[ligname][2]:
                dent_i = 1
            else:
                if isinstance(licores[ligname][2], (str, unicode)):
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
    # sort by descending denticity (needed for adjacent connection atoms)
    ligandsU,occsU,dentsU = ligs,occs0,dentl # save unordered lists
    indcs = smartreorderligs(args,ligs,dentl,licores)
    ligands = [ligs[i] for i in indcs]  # sort ligands list
    occs = [occs0[i] for i in indcs]    # sort occurrences list
    tcats = [cats0[i] for i in indcs]   # sort connections list
    dents = [dentl[i] for i in indcs]   # sort denticities list
    # if using decorations, make repeatable list
    if args.decoration:
        if not args.decoration_index:
            print('Warning, no deocoration index given, assuming first ligand')
            args.decoration_index = [[0]]
        if len(args.decoration_index) != len(ligs):
            new_decoration_index =  []
            new_decorations = []
            for i in range(0,len(ligs)):
                if len(args.decoration_index) > i:
                    new_decoration_index.append(args.decoration_index[i])
                    new_decorations.append(args.decoration[i])
                else:
                    new_decoration_index.append([])
                    new_decorations.append(False)
            if args.debug:
                print('setting decoration:')
                print(new_decoration_index)
                print(new_decorations)
            args.decoration = new_decorations
            args.decoration_index =  new_decoration_index
        args.decoration_index = [args.decoration_index[i] for i in indcs]   # sort decorations list
        args.decoration = [args.decoration[i] for i in indcs]   # sort decorations list
    # sort keepHs list and unpack into list of tuples representing each connecting atom###
    keepHs = [k for k in args.keepHs]
    keepHs = [keepHs[i] for i in indcs]
    for i,keepH in enumerate(keepHs):
        keepHs[i] = [keepHs[i]] * dents[i]
    # sort M-L bond list
    MLb = False
    if args.MLbonds:
        MLb = [k for k in args.MLbonds]
        for j in range(len(args.MLbonds),len(ligs)):
            MLb.append(False)
        MLb = [MLb[i] for i in indcs] # sort MLbonds list
    # sort ligands custom angles
    pangles = False
    if args.pangles:
        pangles = []
        for j in range(len(args.pangles),len(ligs)):
            pangles.append(False)
        pangles = [args.pangles[i] for i in indcs] # sort custom langles list

    # compute number of connecting points required
    cpoints_required = 0
    for i,ligand in enumerate(ligands):
        for j in range(0,occs[i]):
            cpoints_required += dents[i]

    # load core and initialize template
    m3D,core3D,geom,backbatoms,coord,corerefatoms = init_template(args,cpoints_required,globs)
    # Get connection points for all the ligands
    # smart alignment and forced order

    #if geom:
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
    if not geom:
        for comb in batslist:
            for i in comb:
                if i == 1:
                    batslist[comb][i] = m3D.natoms - coord + 1
    # initialize ANN
    ANN_flag,ANN_bondl,ANN_reason,ANN_attributes = init_ANN(args,ligands,occs,dents,batslist,tcats,licores)
    this_diag.set_ANN(ANN_flag,ANN_reason,ANN_attributes)

    # freeze core
    for i in range(0,core3D.natoms):
        frozenats.append(i)

    # loop over ligands and begin functionalization
    # loop over ligands
    totlig = 0  # total number of ligands added
    ligsused = 0
    for i,ligand in enumerate(ligands):
        if not(ligand=='x' or ligand =='X'):
            # load ligand
            lig,emsg = lig_load(ligand)
            # add decorations to ligand
            if args.decoration and args.decoration_index:
                if len(args.decoration) > i and len(args.decoration_index) > i:
                    if args.decoration[i]:
                        if args.debug:
                            print('decorating ' + str(ligand) + ' with ' +str(args.decoration[i]) + ' at sites '  + str(args.decoration_index))
                        lig = decorate_ligand(args,ligand,args.decoration[i],args.decoration_index[i])
                    else:
                        #keeps ligands that are not being decorated
                        lig,emsg = lig_load(ligand)
            lig.convert2mol3D()
            # initialize ligand
            lig3D,rempi,ligpiatoms = init_ligand(args,lig,tcats,keepHs,i)
            if emsg:
                return False,emsg
        for j in range(0,occs[i]):
            denticity = dents[i]
            if not(ligand=='x' or ligand =='X') and (totlig-1+denticity < coord):
                # add atoms to connected atoms list
                catoms = lig.cat # connection atoms
                initatoms = core3D.natoms # initial number of atoms in core3D
                for at in catoms:
                    connected.append(initatoms+at)
                # initialize variables
                mcoords = core3D.getAtom(0).coords() # metal coordinates in backbone
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0 # initialize variables
                coreref = corerefatoms.getAtom(totlig)
                # connecting point in backbone to align ligand to
                batoms = get_batoms(args,batslist,ligsused)
                cpoint = m3D.getAtom(batoms[0])
                # attach ligand depending on the denticity
                # optimize geometry by minimizing steric effects
                print batoms
                if (denticity == 1):
                    lig3D,MLoptbds = align_dent1_lig(args,cpoint,core3D,coreref,ligand,lig3D,catoms,rempi,ligpiatoms,MLb,ANN_flag,ANN_bondl,this_diag,MLbonds,MLoptbds,i)
                elif (denticity == 2):
                    lig3D,frozenats,MLoptbds = align_dent2_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl,this_diag,MLbonds,MLoptbds,frozenats,i)
                elif (denticity == 3):
                    lig3D,frozenats,MLoptbds = align_dent3_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl,this_diag,MLbonds,MLoptbds,frozenats,i)
                elif (denticity == 4):
					# note: catoms for ligand should be specified clockwise
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
                    # align ligand center of symmetry to core reference atom
                    auxmol_lig = mol3D()
                    auxmol_m3D = mol3D()
                    for iiax in range(0,4):
                        auxmol_lig.addAtom(lig3D.getAtom(catoms[iiax]))
                        auxmol_m3D.addAtom(m3D.getAtom(batoms[iiax]))
                    lig3D.alignmol(atom3D('C',auxmol_lig.centersym()),m3D.getAtom(0))
                    # necessary to prevent lig3D from being overwritten
                    lig3Dtmp = mol3D()
                    lig3Dtmp.copymol3D(lig3D)
                    # compute average metal-ligand distance
                    auxmol_lig = mol3D()
                    auxmol_m3D = mol3D()
                    sum_MLdists = 0
                    for iiax in range(0,4):
                        auxmol_lig.addAtom(lig3Dtmp.getAtom(catoms[iiax]))
                        auxmol_m3D.addAtom(m3D.getAtom(batoms[iiax]))
                        sum_MLdists += distance(m3D.getAtomCoords(0),auxmol_lig.getAtomCoords(iiax))
                    avg_MLdists = sum_MLdists/4
                    # scale template by average M-L distance
                    auxmol_m3D.addAtom(m3D.getAtom(0))
                    for iiax in range(0,4):
                        auxmol_m3D.BCM(iiax,4,avg_MLdists)
                    auxmol_m3D.deleteatom(4)
                    # align lig3D to minimize RMSD from template
                    auxmol_lig,U,d0,d1 = kabsch(auxmol_lig,auxmol_m3D)
                    lig3D.translate(d0)
                    lig3D = rotate_mat(lig3D,U)
                    ## align plane
                    #r0c = m3D.getAtom(batoms[0]).coords()
                    #r1c = m3D.getAtom(batoms[1]).coords()
                    #r2c = m3D.getAtom(batoms[2]).coords()
                    #r0l = lig3D.getAtom(catoms[0]).coords()
                    #r1l = lig3D.getAtom(catoms[1]).coords()
                    #r2l = lig3D.getAtom(catoms[2]).coords()
                    #theta,uc = rotation_params(r0c,r1c,r2c) # normal vector to backbone plane
                    #theta,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    #lig3Db = mol3D()
                    #lig3Db.copymol3D(lig3D)
                    #theta = 180*arccos(dot(uc,ul)/(norm(uc)*norm(ul)))/pi
                    #u = cross(uc,ul)
                    ## rotate around axis to match planes
                    #theta = 180-theta if theta > 90 else theta
                    #lig3D = rotate_around_axis(lig3D,r0l,u,theta)

                    ## rotate ar?ound secondary axis to match atoms
                    #r0l = lig3D.getAtom(catoms[0]).coords()
                    #r1l = lig3D.getAtom(catoms[1]).coords()
                    #r2l = lig3D.getAtom(catoms[2]).coords()
                    #theta0,ul = rotation_params(r0l,r1l,r2l) # normal vector to ligand plane
                    #rm = lig3D.centersym()
                    #r1 = vecdiff(r0l,mcoords)
                    #r2 = vecdiff(r0c,mcoords)
                    #theta = 180*arccos(dot(r1,r2)/(norm(r1)*norm(r2)))/pi
                    #lig3Db = mol3D()
                    #lig3Db.copymol3D(lig3D)
                    #if args.debug:
                        #print('normal to tetradentate ligand plane: ',ul)
                        #print('lig center of symm ',rm)
                        #lig3D.writexyz('lig3d.xyz')
                        #lig3Db.writexyz('lig3db.xyz')
                    ## rotate around axis and get both images
                    #lig3D = rotate_around_axis(lig3D,mcoords,ul,theta)

                    bondl = get_MLdist(args,lig3D,atom0,ligand,m3D.getAtom(0),MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
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
                    lig3D.translate(vecdiff(mcoords,ligc.centersym()))
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
                    bondl,exact_match = get_MLdist_database(args,core3D.getAtom(0),lig3D,catoms[0],ligand,MLbonds)
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
                    core3D.translate(vecdiff(ligc.centersym(),mcoords))
                    bondl,exact_match = get_MLdist_database(args,core3D.getAtom(0),lig3D,catoms[0],ligand,MLbonds)
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
                core3D.convert2OBMol()
                core3D.convert2mol3D()
                # remove dummy cm atom if requested
                if rempi:
                    core3D.deleteatom(core3D.natoms-1)
                if args.calccharge:
                    core3D.charge += lig3D.charge
                # perform FF optimization if requested

                if 'a' in args.ffoption:
                    print('FF optimizing molecule after placing ligand')
                    core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
            totlig += denticity
            ligsused += 1
    # perform FF optimization if requested
    if 'a' in args.ffoption:
        print('Performing final FF opt')
        core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
    return core3D,complex3D,emsg,this_diag


## Main structure generation routine - single structure
#  @param strfiles List of xyz files generated
#  @param args Namespace of arguments
#  @param rootdir Directory of current run
#  @param ligands List of ligands
#  @param ligoc List of ligand occupations
#  @param globs Global variables
#  @param sernum Serial number of complex for naming
#  @param nconf Conformer ID, if any
#  @return List of xyz files generated, error messages
def structgen_one(strfiles,args,rootdir,ligands,ligoc,globs,sernum,nconf=False):
    # load ligand dictionary
    licores = getlicores()
    # build structure
    sanity = False
    this_diag = run_diag()
    if (ligands):
        core3D,complex3D,emsg,this_diag = mcomplex(args,ligands,ligoc,licores,globs)
        name_core = args.core
        if emsg:
            return False,emsg
    else:
        print('You specified no ligands. Returning the core.')
        core3D = mol3D()
        name_core = core3D
    # generate file name parts
    ligname = ''
    nosmiles = 0
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
    if args.calccharge:
        args.charge = core3D.charge
        if args.debug:
            print('setting charge to be ' + str(args.charge))
    # check for molecule sanity
    sanity,d0 = core3D.sanitycheck(True)
    if sanity:
        print 'WARNING: Generated complex is not good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
        if args.gui:
            ssmsg = 'Generated complex in folder '+rootdir+' is no good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
            qqb = mQDialogWarn('Warning',ssmsg)
            qqb.setParent(args.gui.wmain)
    if args.debug:
        print('setting sanity diag, min dist at ' +str(d0) + ' (higher is better)')
    this_diag.set_sanity(sanity,d0)
    # generate file name
    fname = name_complex(rootdir,name_core,ligands,ligoc,sernum,args,nconf,sanity)
    # write xyz file
    core3D.writexyz(fname)
    strfiles.append(fname)
    # write report file
    this_diag.set_mol(core3D)
    this_diag.write_report(fname+'.report')
    # write input file from command line arguments
    getinputargs(args,fname)
    del core3D
    return strfiles, emsg, this_diag

## Main structure generation routine - multiple structures
#  @param args Namespace of arguments
#  @param rootdir Directory of current run
#  @param ligands List of ligands
#  @param ligoc List of ligand occupations
#  @param globs Global variables
#  @param sernum Serial number of complex for naming
#  @return List of xyz files generated, error messages
def structgen(args,rootdir,ligands,ligoc,globs,sernum):
    emsg = False
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn

    strfiles = []

    if args.smicat:
        if sum([len(i)>1 for i in args.smicat]) > 0:
            print('You have specified multidentate SMILES ligand(s).')
            print('We will automatically find suitable conformer(s) for coordination.')
            for n in range(1,int(args.nconfs)+1):
                print 'Generating conformer '+str(n)+' of '+args.nconfs+':'
                strfiles, emsg, this_diag = structgen_one(strfiles,args,rootdir,ligands,ligoc,globs,sernum,n)
                print strfiles
        else:
            strfiles, emsg, this_diag = structgen_one(strfiles,args,rootdir,ligands,ligoc,globs,sernum)
    else:
        strfiles, emsg, this_diag = structgen_one(strfiles,args,rootdir,ligands,ligoc,globs,sernum)

    # score conformers
    conf3Ds = dict()
    if args.smicat:
        if args.scoreconfs:
            for i,strfile in enumerate(strfiles):
                fname = strfile.rsplit('/',-1)[-1]
                print rootdir
                conf3D = core_load(strfile+'/'+fname+'.xyz')
                conf3Ds[i] = conf3D
    #print conf3Ds

    # number of different combinations
    if args.bindnum and args.bind:
        Nogeom = int(args.bindnum)
    elif args.smicat:
        if sum([len(i)>=1 for i in args.smicat]) > 0:
            Nogeom = int(args.nconfs)
    else:
        Nogeom = 1
    # generate multiple geometric arrangements
    if args.bind:
        # load bind, add hydrogens and convert to mol3D
        bind,bsmi,emsg = bind_load(args.bind)
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
        # fetch base name
        fname = get_name(args,rootdir,core,ligname,bind,bsmi)
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
                    if ('a' in args.place):
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

    pfold = rootdir.split('/',1)[-1]
    if args.gui:
        args.gui.iWtxt.setText('In folder '+pfold+' generated '+str(Nogeom)+' structures!\n'+args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    print '\nIn folder '+pfold+' generated ',Nogeom,' structure(s)!'

    return strfiles, emsg, this_diag
