## file structgen.py
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
from molSimplify.Scripts.krr_prep import *
import os, sys, time
from pkg_resources import resource_filename, Requirement
import openbabel, random, itertools, numpy
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
    globs = globalvars()
    catalysis_flag = False
    if args.skipANN:
         print('Skipping ANN')
         ANN_flag = False
         ANN_bondl = len([item for items in batslist for item in items])*[False] ## there needs to be 1 length per possible lig
         ANN_reason = 'ANN skipped by user'
    else:

         try:
         # if True:
            if args.oldANN:
                print('using old ANN by request')
                ANN_flag,ANN_reason,ANN_attributes = ANN_preproc(args,ligands,occs,dents,batslist,tcats,licores)
            else:
                if globs.testTF():
                    ## new RACs-ANN
                    from molSimplify.Scripts.tf_nn_prep import tf_ANN_preproc
                    ANN_flag,ANN_reason,ANN_attributes, catalysis_flag = tf_ANN_preproc(args,ligands,occs,dents,batslist,tcats,licores)
                else:
                    # old MCDL-25
                    print('using old ANN because tensorflow/keras import failed')
                    ANN_flag,ANN_reason,ANN_attributes = ANN_preproc(args,ligands,occs,dents,batslist,tcats,licores)
            if ANN_flag:
                 ANN_bondl = ANN_attributes['ANN_bondl']
                 print('ANN bond length is ' + str(ANN_bondl))
            else:
                 ANN_bondl = len([item for items in batslist for item in items])*[False] ## there needs to be 1 length per possible lig
                 if args.debug:
                     print("ANN called failed with reason: " + ANN_reason)
         except:
         # else:
             print("ANN call rejected")
             ANN_reason = 'uncaught exception'
             ANN_flag = False
             ANN_bondl =  len([item for items in batslist for item in items])*[False]
    return ANN_flag,ANN_bondl,ANN_reason,ANN_attributes, catalysis_flag


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

## Initializes core and template mol3Ds and properties
#  @param args Namespace of arguments
#  @param cpoints_required Number of connecting points required
#  @return mol3D of core, template, geometry, backbone atoms, coordination number, core reference atom index
def init_mcomplex_template(args, core3D, cpoints_required, mligcatoms_ext, bondlss, bangless, bdihedralss):
    # initialize core and template
    m3D = mol3D()
    # container for ordered list of core reference atoms
    corerefatoms = mol3D()
    # geometry load flag
    geom = False
    cpoint_nums = []
    backbatoms = []
    # bring mligcatoms_ext that are corresponding to core to the front
    i = 0
    # print(mligcatoms_ext)
    for idx, mligcatom_ext in enumerate(mligcatoms_ext):
        mligcatomsym = core3D.getAtom(mligcatom_ext).sym
        if mligcatomsym in [core.lower() for core in args.core]:
            mligcatoms_ext[idx], mligcatoms_ext[i] = mligcatoms_ext[i], mligcatoms_ext[idx]
            bondlss[idx], bondlss[i] = bondlss[i], bondlss[idx]
            bangless[idx], bangless[i] = bangless[i], bangless[idx]
            i += 1
    # check mligcatoms
    d_ref = 100
    m3D.copymol3D(core3D)
    for i in range(cpoints_required):
        bondl_core3D, bondl_m3D, bondl_sub = bondlss[i][0], bondlss[i][1], bondlss[i][2]
        bangle_m3D, bangle_core3D, bangle_sub = bangless[i][0], bangless[i][2], bangless[i][3]
        bdihedral_m3D, bdihedral_core3D, bdihedrel_sub = bdihedralss[i][0], bdihedralss[i][2], bdihedralss[i][3]
        mligcatom_ext = mligcatoms_ext[i]
        mligcatomcoords = m3D.getAtom(mligcatom_ext).coords()
        mligcatomsym = m3D.getAtom(mligcatom_ext).sym
        midx = m3D.findMetal()[0]
        mcoord = m3D.getAtom(midx).coords()
        ## adjust mlig
        # obtain the anchoring atom
        for idx in m3D.getBondedAtoms(mligcatom_ext):
            coord = m3D.getAtom(idx).coords()
            d = distance(mcoord, coord)
            if d < d_ref:
                d_ref = d
                mligcatom_ext_anchor = idx
                mliganchorcoords = coord
        # adjust M-L distance for the incoming substrate
        # if (args.core not in m3D.getAtom(mligcatom_ext).sym):
        core3D.BCM(mligcatom_ext, mligcatom_ext_anchor, bondl_core3D)
        m3D.BCM(mligcatom_ext, mligcatom_ext_anchor, bondl_core3D)
        # check natoms in mlig
        if args.mlig:
            lig, emsg = lig_load(args.mlig[0])
            lig.convert2mol3D()
            natoms_lig = lig.natoms
            # adjust M-{L} angle for the incoming substrate
            if natoms_lig > 1:
                core3D.ACM(mligcatom_ext, mligcatom_ext_anchor, midx, bangle_core3D)
                m3D.ACM(mligcatom_ext, mligcatom_ext_anchor, midx, bangle_core3D)
        # if mligcatomsym in [core.lower() for core in args.core]:
        #
        if i == 0:
            # print('mligcatom_ext_anchor is ' + str(mligcatom_ext_anchor))
            # print('bangle_m3D is ' + str(bangle_m3D))
            # cpoint = getconnectiongivenphi(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D)
            cpoint = getconnectiongivenangles(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, (bangle_m3D, bdihedral_m3D))
            # cpoint = getconnection(m3D, mligcatom_ext, bondl_m3D)
            # store core reference atom
            conatom3D = atom3D(m3D.getAtom(mligcatom_ext).sym, m3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print(corerefatoms.getAtom(0).symbol())
            # initiate dummy atom
            dummy_atom = atom3D(Sym='X',xyz=cpoint)
            m3D.addAtom(dummy_atom)
            cpoint_nums.append(m3D.natoms - 1)
            # nums = m3D.findAtomsbySymbol('I')
            backbatoms = getbackbcombsall(cpoint_nums)
            mligcatom_old = mligcatom_ext
        # if the second cpoint is on the same mligcatom as the first one
        elif mligcatom_ext == mligcatom_old:
            cpoint = getconnectiongivenr(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D, bondl_sub)
            # print('bondl_sub is ' + str(bondl_sub))
            # store core reference atom
            conatom3D = atom3D(core3D.getAtom(mligcatom_ext).sym, core3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print(corerefatoms.getAtom(0).symbol())
            #corerefatoms.append(ccatoms[i])N
            # add connecting points to template
            m3D.addAtom(atom3D(Sym='X', xyz=cpoint))
            cpoint_nums.append(m3D.natoms - 1)
            # except IndexError:
            #     pass
            backbatoms = getbackbcombsall(cpoint_nums)
        else:
            cpoint = getconnectiongivenangles(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, (bangle_m3D, bdihedral_m3D))
            # cpoint = getconnectiongivenphi(m3D, mligcatom_ext, mligcatom_ext_anchor, bondl_m3D, bangle_m3D)
            # store core reference atom
            conatom3D = atom3D(core3D.getAtom(mligcatom_ext).sym, core3D.getAtom(mligcatom_ext).coords())
            corerefatoms.addAtom(conatom3D)
            if args.debug:
                print(corerefatoms.getAtom(0).symbol())
            # initiate dummy atom
            dummy_atom = atom3D(Sym='X', xyz=cpoint)
            m3D.addAtom(dummy_atom)
            cpoint_nums.append(m3D.natoms - 1)
            # obtain rotation axis
            mliganchorcoords = m3D.getAtomCoords(mligcatom_ext_anchor)
            mligcatoms_ext_anchor2 = [idx for idx in m3D.getBondedAtoms(mligcatom_ext_anchor) if idx != mligcatom_ext]
            for idx in mligcatoms_ext_anchor2:
                coord = m3D.getAtom(idx).coords()
                d = distance(mcoord, coord)
                mliganchor2coords = coord
                if d < d_ref:
                    d_ref = d
                    mliganchor2coords = coord
            r_mliganchor_mliganchor2 = vecdiff(mliganchorcoords, mliganchor2coords)
            # rotate the [mligcatom_ext, mligcatom_ext_anchor] bond to bring the dummy atoms closer
            d0 = 100
            refidx = sorted(m3D.findAtomsbySymbol('X'))[0]
            refcoords = m3D.getAtomCoords(refidx)
            num_rotation = 0
            m3D_ = mol3D()
            m3D_.copymol3D(m3D)
            theta = 5
            while num_rotation < 72:
                num_rotation += 1
                m3D_.ACM_axis(mligcatom_ext, mligcatom_ext_anchor, r_mliganchor_mliganchor2, theta)
                xidx = sorted(m3D_.findAtomsbySymbol('X'))[-1]
                xcoords = m3D_.getAtomCoords(xidx)
                d = distance(refcoords, xcoords)
                if d < d0:
                    d0 = d
                    m3D = mol3D()
                    m3D.copymol3D(m3D_)
    # set charge from oxidation state if desired
    if args.calccharge:
        if args.oxstate:
            if args.oxstate in romans.keys():
                core3D.charge = int(romans[args.oxstate])
            else:
                core3D.charge = int(args.oxstate)
    # remove X atoms from m3D to generate core3D
    core3D = mol3D()
    core3D.copymol3D(m3D)
    for atidx in range(core3D.natoms)[::-1]:
        asym = core3D.getAtom(atidx).sym
        if asym == 'X':
            core3D.deleteatom(atidx)

    return m3D, core3D, geom, backbatoms, coord, corerefatoms

## Initializes substrate 3D geometry and properties
#  @param args Namespace of arguments
#  @param sub mol3D of ligand
#  @param subcatoms list of connecting atoms in the substrate
#  @return mol3D of ligand, flag for pi-coordination, pi-coordinating atoms
def init_substrate(args, sub, subcatoms, bondlss, bangless, bdihedralss):
    rempi = False
    idx3 = ''
    sub.convert2mol3D()
    sub3D = mol3D()
    sub3D.copymol3D(sub)
    # sub3D.createMolecularGraph(False)
    for rxn_type_i in range(len(bondlss)):
        bondl_m3D = bondlss[rxn_type_i][1]
        bondl_sub = bondlss[rxn_type_i][2]
        bangle_sub = bangless[rxn_type_i][3]
        bdihedral_sub = bdihedralss[rxn_type_i][3]
        # if SMILES string, copy connecting atoms list to mol3D properties
        # if not sub.cat and tcats[i]:
        #     if 'c' in tcats[i]:
        #         lig.cat = [lig.natoms]
        #     else:
        #         lig.cat = tcats[i]
        # change name
        # substrate bond manipulation to prepare for activation
        anchor_atom_idx = int(subcatoms[rxn_type_i])
        ref_to_sub = False
        moved_atom_idx = ''
        for atidx in sub3D.getBondedAtoms(anchor_atom_idx):
            moved_atom_coords = sub3D.getAtomCoords(atidx)
            # check to see if there is another connection atom
            distss = []
            if atidx in subcatoms:
                moved_atom_idx = atidx
                break
            else:
                ref_to_sub = True
                for subcatom in subcatoms:
                    subcoords = sub3D.getAtomCoords(subcatom)
                    d = distance(moved_atom_coords, subcoords)
                    distss.append((d, atidx))
        if ref_to_sub:
            distss.sort()
            moved_atom_idx = distss[0][1]
        elif not moved_atom_idx:
            moved_atom_idx = sub3D.getBondedAtoms(anchor_atom_idx)[0]
        sub3D.BCM(moved_atom_idx, anchor_atom_idx, bondl_sub)
        # angle adjustment
        if sub3D.natoms > 2 and len(sub3D.getBondedAtoms(moved_atom_idx)) > 1 and moved_atom_idx not in subcatoms:
            idx3 = [idx for idx in sub3D.getBondedAtoms(moved_atom_idx) if idx != anchor_atom_idx][0]
            sub3D.ACM(anchor_atom_idx, moved_atom_idx, idx3, bangle_sub)
        # # check for pi-coordinating substrate
        subpiatoms = []
        if 'pi' in sub.cat:
            if args.debug:
                print('substrate is identified as a pi-type substrate')
            subcatoms = []
            for k in sub.cat[:-1]:
                if isinstance(k,int):
                    subcatoms.append(k)
                # sub3Dpiatoms.addAtom(sub3D.getAtom(k))
                # sub3Dpiatoms.addAtom(sub3D.getAtom(k))
            # subpiatoms = sub.cat[:-1]
            # sub3D.addAtom(atom3D('C',sub3Dpiatoms.centermass()))
            # if args.debug:
            #     sub3D.printxyz()
            # sub3D.cat = [sub3D.natoms-1]
            rempi = True
        cpoint = getconnectiongivenangles(sub3D, anchor_atom_idx, moved_atom_idx, bondl_m3D, (bangle_sub, bdihedral_sub))
        dummy = atom3D(Sym='X', xyz=cpoint)
        sub3D.addAtom(dummy)
        xidxes = sorted(sub3D.findAtomsbySymbol('X'))
        refidx = xidxes[0]
        refcoords = sub3D.getAtomCoords(refidx)
        if len(xidxes) > 1 and idx3:
            anchor_coords = sub3D.getAtomCoords(anchor_atom_idx)
            moved_coords = sub3D.getAtomCoords(moved_atom_idx)
            idx3_coords = sub3D.getAtomCoords(idx3)
            u0 = vecdiff(anchor_coords, moved_coords)
            u1 = vecdiff(moved_coords, idx3_coords)
            # u =
            d0 = 10
            thetas = [0, 0]
            theta1 = 5
            for theta0 in range(0, 360, 5):
                num_rotation = 0
                sub3D_ = mol3D()
                sub3D_.copymol3D(sub3D)
                sub3D_.ACM_axis(anchor_atom_idx, moved_atom_idx, u0, theta0)
                while num_rotation < 72:
                    num_rotation += 1
                    sub3D_.ACM_axis(anchor_atom_idx, moved_atom_idx, u1, theta1)
                    pidx = xidxes[1]
                    pcoords = sub3D_.getAtomCoords(pidx)
                    d = distance(refcoords, pcoords)
                    # planar = checkplanar(anchor_coords, moved_coords, pcoords, refcoords)
                    if d < d0:
                        d0 = d
                        thetas[0] = theta0
                        thetas[1] = theta1 * num_rotation
            # print('distance is ' + str(d0))
            # print('theta is ' + str(thetas))
            sub3D.ACM_axis(anchor_atom_idx, moved_atom_idx, u0, thetas[0])
            sub3D.ACM_axis(anchor_atom_idx, moved_atom_idx, u1, thetas[1])

    return sub3D, rempi, subcatoms, subpiatoms

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
            lig3D.convert2mol3D()
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

    if lig.needsconformer:
        tcats[i] = True
        print('getting conformers for ' + str(lig.ident))

    if len(lig.cat) > 1 and tcats[i]:
        print('generating conformations')
        # loop  over conformation gen until success or break 
        breaker = False
        count = 0
        while (not breaker) and count <= 5:
            try:
                lig3D = GetConf(lig3D,lig.cat)    # check if ligand should decorated
                breaker = True
            except:
                count += 1
                print('lig conformer input failed ' + str(count)  + ' times, trying again...')
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
        ff = 'uff'
    if debug:
        print('using ff: ' + ff)
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
                if debug:
                    print('using connnected opt to freeze atom number: ' + str(catom))
            else:
                constr.AddDistanceConstraint(midx+1,catom+1,mlbonds[ii]) # indexing babel
        #print('ff is '+ str(ff))
        if not ff.lower() == "uff":
            bridgingatoms = []
            # identify bridging atoms in the case of bimetallic cores, as well as single-atom ligands (oxo, nitrido)
            # these are immune to deletion
            for i in range(mol.natoms):
                nbondedmetals = len([idx for idx in range(len(mol.getBondedAtoms(i))) if mol.getAtom(mol.getBondedAtoms(i)[idx]).ismetal()])
                if nbondedmetals > 1 or (nbondedmetals == 1 and len(mol.getBondedAtoms(i)) == 1):
                    bridgingatoms.append(i)
            # ensure correct valences for FF setup
            deleted_bonds = 0

            for m in indmtls:
            # first delete all metal-ligand bonds excluding bridging atoms
                for i in range(len(mol.getBondedAtoms(m))):
                    if OBMol.GetBond(m+1,mol.getBondedAtoms(m)[i]+1) is not None and mol.getBondedAtoms(m)[i] not in bridgingatoms:
                        OBMol.DeleteBond(OBMol.GetBond(m+1,mol.getBondedAtoms(m)[i]+1))
                        #print('FFopt deleting bond')
                        deleted_bonds += 1
                print('FFopt deleted ' +str(deleted_bonds) + ' bonds')
                # then add back one metal-ligand bond for FF
                if OBMol.GetAtom(m+1).GetValence() == 0:
                    for i in mol.getBondedAtoms(m):#getBondedAtomsOct(m,deleted_bonds+len(bridgingatoms)):
                        if OBMol.GetAtom(m+1).GetValence() < 1 and i not in bridgingatoms:
                            OBMol.AddBond(m+1,i+1,1)
        # freeze small ligands
        for cat in frozenats:
            if debug:
                print('using frozenats to freeze atom number: ' + str(cat))
            constr.AddAtomConstraint(cat+1) # indexing babel
        if debug:
            for iiat,atom in enumerate(openbabel.OBMolAtomIter(OBMol)):
                print (' atom '+str(iiat)+' atomic num '+str(atom.GetAtomicNum())+' valence '+str(atom.GetValence())+ ' is fixed '+ str(constr.IsFixed(iiat+1)))
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
            if debug:
                print('running ' +str(n) + ' steps')
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
            print('doing 50 steps')
            forcefield.ConjugateGradients(50)
        else:
            print('doing 200 steps')
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

## Finds the optimum attachment point for an atom/group to a central atom given the desired bond length and bond angle
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param refidx idx of the atom bonded to the core connecting atom
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnectiongivenphi(core,cidx,refidx,BL,BA):
    # ncore = core.natoms
    # groups = core.getBondedAtoms(cidx)
    ccoords = core.getAtom(cidx).coords()
    refcoords = core.getAtomCoords(refidx)
    r_c_ref = vecdiff(refcoords, ccoords)
    pcoords = getPointu(ccoords, BL, r_c_ref)
    # print('BL is ' + str(BL))
    theta, u = rotation_params(refcoords, ccoords, pcoords)
    phi = BA / 180. * pi
    # pcoords = PointTranslateSphgivenphi(ccoords, pcoords,[BL, phi])
    # pcoords = PointRotateAxis(u, ccoords, pcoords, phi)
    # print('pc distance is ' + str(distance(pcoords, ccoords)))
    # print(pcoords)
    # if core.getAtom(cidx).sym in [core.getAtom(midx).sym]:
    #     refcoords = [0,0,0]
    #     for fidx in core.getBondedAtoms(cidx):
    #         nfidx = len(core.getBondedAtoms(cidx))
    #         refcoords[0] += core.getAtom(fidx).coords()[0]/nfidx
    #         refcoords[1] += core.getAtom(fidx).coords()[1]/nfidx
    #         refcoords[2] += core.getAtom(fidx).coords()[2]/nfidx
    # else:
    #     refcoords = core.getAtom(refidx).coords()
    # phi = float(BA/180.*pi)
    # brute force search
    cpoint = []
    objopt = 0
    for itheta in range(1,359,1):
        theta = itheta / 180. * pi
        pcoords_ = PointTranslateSph(ccoords, pcoords,[BL, theta, phi])
        # P = PointTranslateSph(ccoords,ccoords,[BL,itheta,BA])
        dists = []
        for ig in range(core.natoms):
            if ig != cidx:
                dists.append(distance(core.getAtomCoords(ig), pcoords_))
        obj = min(dists)
        if obj > objopt:
            # print('obj is ' + str(obj))
            objopt = obj
            cpoint = pcoords_
    # P = PointTranslateSphgivenphi(ccoords,refcoords,[BL,phi])
    # cpoint = pcoords_

    return cpoint

## Finds the optimum attachment point for an atom/group to a central atom given the desired bond length and bond angle
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param refidx idx of the atom bonded to the core connecting atom
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnectiongivenangles(core, cidx, refidx, BL, BAs):
    # ncore = core.natoms
    # groups = core.getBondedAtoms(cidx)
    cpoint = []
    ccoords = core.getAtom(cidx).coords()
    refcoords = core.getAtomCoords(refidx)
    anchoridxes = [idx for idx in core.getBondedAtoms(cidx) if idx != refidx]
    # print('cidx is ' + str(cidx))
    # print('refidx is ' + str(refidx))
    # print('anchoridxes are ' + str(anchoridxes))
    r_c_ref = vecdiff(refcoords, ccoords)
    pcoords = getPointu(ccoords, BL, r_c_ref)
    # print('BL is ' + str(BL))
    # theta, u = rotation_params(refcoords, ccoords, pcoords)
    phi = BAs[0] / 180. * pi
    phi_ref = BAs[0]
    theta_ref = BAs[1]
    objopt = 0
    if anchoridxes:
        # print('theta_ref is ' + str(theta_ref))
        anchorcoords = core.getAtomCoords(anchoridxes[0])
        theta, norm = rotation_params(refcoords, ccoords, anchorcoords)
        # brute force search
        objopt = 0
        distss = []
        # print('phi_ref is ' + str(phi_ref))
        # print('theta_ref is ' + str(theta_ref))
        #! TODO this double for loop needs to go
        for itheta in range(0, 360, 10):
            theta = itheta / 180. * pi
            for iphi in range(0, 360, 10):
                phi = iphi / 180. * pi
                pcoords_ = PointTranslateSph(ccoords, r_c_ref, [BL, theta, phi])
                angle = vecangle(norm, r_c_ref)
                angle_phi = vecangle(r_c_ref, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                if angle < 80:
                    angle_theta = theta_ref
                else:
                    angle_theta = vecangle(norm, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                # print('angle is ' + str(angle))
                if abs(angle_theta - theta_ref) < 5 and abs(angle_phi - phi_ref) < 5:
                    dists = []
                    xidxes = core.findAtomsbySymbol('X')
                    # whether the second X should be placed near the first one
                    if xidxes:
                        for xidx in xidxes:
                            dists.append(-1 * distance(core.getAtomCoords(xidx), pcoords_))
                    else:
                        for ig in range(core.natoms):
                            if ig != cidx:
                                dists.append(distance(core.getAtomCoords(ig), pcoords_))
                    distss.append((min(dists), pcoords_))
                    # print('distss are ' + str(distss))
                    # print('angle_phi is ' + str(angle_phi) + ', angle_theta is ' + str(angle_theta))
                core3D = mol3D()
                core3D.copymol3D(core)
                atom = atom3D(Sym='X', xyz=pcoords_)
                core3D.addAtom(atom)
        # print('distss are ' + str(sorted(distss)[-1]))
        cpoint = sorted(distss)[-1][1]
    else:
        # print('theta_ref is ' + str(theta_ref))
        distss = []
        for itheta in range(1, 359, 1):
            theta = itheta / 180. * pi
            for iphi in range(0, 360, 10):
                phi = iphi / 180. * pi
                pcoords_ = PointTranslateSph(ccoords, pcoords, [BL, theta, phi])
                angle_phi = vecangle(r_c_ref, [pcoords_[i] - ccoords[i] for i in range(len(ccoords))])
                if abs(angle_phi - phi_ref) < 5:
                    dists = []
                    for ig in range(core.natoms):
                        if ig != cidx:
                            dists.append(distance(core.getAtomCoords(ig), pcoords_))
                    distss.append((min(dists), pcoords_))
        cpoint = sorted(distss)[-1][1]

    return cpoint

## Finds the optimum attachment point for an atom/group to a central atom given the desired bond length and bond angle
#
#  Objective function maximizes the minimum distance between attachment point and other groups bonded to the central atom
#  @param core mol3D of core
#  @param cidx Core connecting atom index
#  @param refidx idx of the atom bonded to the core connecting atom
#  @param BL Optimal core-ligand bond length
#  @return Coordinates of optimum attachment point
def getconnectiongivenr(core,cidx,refidx,BL,BA,r):
    midx = core.findMetal()[0]
    ccoords = core.getAtom(cidx).coords()
    if core.getAtom(cidx).sym in [core.getAtom(midx).sym]:
        refcoords = [0,0,0]
        for fidx in core.getBondedAtoms(cidx):
            nfidx = len(core.getBondedAtoms(cidx))
            refcoords[0] += core.getAtom(fidx).coords()[0]/nfidx
            refcoords[1] += core.getAtom(fidx).coords()[1]/nfidx
            refcoords[2] += core.getAtom(fidx).coords()[2]/nfidx
    else:
        refcoords = core.getAtom(refidx).coords()
    phi = float(BA/180.*pi)
    # brute force search
    cpoint = []
    P0coords = core.atoms[-1].coords()
    P = PointTranslateSphgivenr(ccoords, refcoords, [BL, phi], P0coords, r)
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
            ## warning: skipping this part because
            ## we no longer understand it
            if False:
                r1 = lig3D.getAtom(atom0).coords()
                r2 = auxmol.getAtom(0).coords()
                theta,u = rotation_params([1,1,1],r1,r2)
                lig3D = rotate_around_axis(lig3D,r1,u,-1*globs.linearbentang)
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

## Aligns a substrate's activated bond along the metal-connecting atom axis
#  @param corerefcoords Core reference coordinates
#  @param sub3D mol3D of substrate
#  @param atom0 Substrate connecting atom index
#  @param core3D mol3D of partially built complex
#  @return mol3D of aligned substrate
def align_sub_firstbond(args,corerefcoords,sub3D,atom0,core3D,bangle_m3Dsub):
    # rotate to align center of symmetry
    # globs = globalvars()
    r0 = corerefcoords
    r1 = sub3D.getAtom(atom0).coords()
    sub3Db = mol3D()
    sub3Db.copymol3D(sub3D)
    # auxmol = mol3D()
    # for at in sub3D.getBondedAtomsSmart(atom0):
    #     auxmol.addAtom(sub3D.getAtom(at))
    at = sub3D.getBondedAtomsSmart(atom0)[0]
    r2 = sub3D.getAtomCoords(at)
    theta,u = rotation_params(r0,r1,r2)
    # rotate around axis and get both images
    sub3D = rotate_around_axis(sub3D,r1,u,bangle_m3Dsub-(180-theta))
    if args.debug:
        print('bangle m3Dsub is ' + str(bangle_m3Dsub))
        print('theta is ' + str(theta))
    sub3D_aligned = mol3D()
    sub3D_aligned.copymol3D(sub3D)

    return sub3D_aligned

## Aligns a linear pi substrate's connecting point to the metal-substrate axis
#  @param corerefcoords Core reference coordinates
#  @param sub3D mol3D of suband
#  @param atom0 substrate connecting atom index
#  @param subpiatoms List of substrate pi-connecting atom indices
#  @return mol3D of aligned substrate
def align_linear_pi_sub(core3D,mligcatoms_ext,sub3D,atom0,subcatoms,bangle_m3D):
    # align the norm of the pi bond to the L-S vector at atom0 (one of the C in the pi bond)
    ratom0 = sub3D.getAtom(atom0).coords()
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords() for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i not in subcatoms][0]
    theta102,uatom102 = rotation_params(ratom1,ratom0,ratom2)
    rL = core3D.getAtom(mligcatoms_ext).coords()
    rLS = vecdiff(rL,ratom0)
    thetaL0u102,uL0u102 = rotation_params(rL,ratom0,uatom102)
    thetaL0u102 = vecangle(rLS,uatom102) # aka u1
    # u2 = cross(uatom102,rLS)
    sub3D_aligned = mol3D()
    sub3D_aligned.copymol3D(sub3D)
    # sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, u2, 90-theta_uatom102_rLS)
    # rotate sub3D such that the substrate lies perpendicular to the complex
    sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, uL0u102, -thetaL0u102)
    # rotate sub3D such that the X-C-M angle matches bangle_m3Dsub
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords() for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i not in subcatoms][0]
    theta10L,u10L = rotation_params(ratom1,ratom0,rL)
    sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, u10L, bangle_m3Dsub+theta10L)
    # agjust the angle among L-C-Ca
    try:
        ratom1 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i in subcatoms][0]
    except:
        ratom1 = [sub3D.getAtom(i).coords() for i in subcatoms if i is not atom0][0]
    ratom2 = [sub3D.getAtom(i).coords() for i in sub3D.getBondedAtoms(atom0) if i not in subcatoms][0]
    # objfuncopt = 90
    # thetaopt = 0
    # for theta in range(0,360,1):
    #     sub3D_aligned = rotate_around_axis(sub3D_aligned,ratom0,rLS, theta)
    #     #objfunc = abs(vecangle(vecdiff(sub3D_aligned.getAtom(atom0).coords(),corerefcoords),vecdiff(sub3D_aligned.getAtom(ligpiatoms[0]).coords(),sub3D_aligned.getAtom(ligpiatoms[1]).coords()))-90)
    #     objfunc = abs(distance(sub3D_aligned.getAtom(ligpiatoms[0]).coords(),corerefcoords) - distance(sub3D_aligned.getAtom(ligpiatoms[1]).coords(),corerefcoords))
    #     if objfunc < objfuncopt:
    #         thetaopt = theta
    #         objfuncopt = objfunc
    #         sub3Dopt = mol3D() # lig3Dopt = sub3D_aligned DOES NOT WORK!!!
    #         sub3Dopt.copymol3D(sub3D_aligned)
    # sub3D_aligned = rotate_around_axis(sub3D_aligned, ratom0, rLS, 30)
    # reduce the steric repulsion between the mcomplex and the substrate by rotating away the ther C in the pi bond
    sub3D_aligned = rotate_MLaxis_minimize_steric(rL,sub3D_aligned,atom0,core3D)

    return sub3D_aligned

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

## Rotates aligned ligand about M-L axis to minimize steric clashes with rest of complex
#  @param mol mol3D of the molecule to be rotated
#  @param molatcoords the coordinates of the atom in the rotated molecule
#  @param refmol mol3D of the molecule of reference
#  @param refcoords the coordinates of the atom in the reference molecule
#  @param u the vector of rotation axis
#  @return mol3D of rotated ligand
def rotate_MLaxis_minimize_steric_ts(mol, coords, refmol, refcoords, u):
    dist0 = 0
    dist = 0
    theta0 = 0
    while theta0 < 360:
        dists = []
        mol_ = rotate_around_axis(mol, refcoords, u, theta0)
        for atom in mol_.atoms:
            coords_ = atom.coords()
            if coords_ != coords or mol_.natoms == 1:
                for atom_ref in refmol.atoms:
                    refcoords_ = atom_ref.coords()
                    if refcoords_ != refcoords:
                        dist = distance(coords_, refcoords_)
                        dists.append(dist)
        dist = min(dists)
        if dist > dist0:
            dist0 = dist
            mol_aligned = mol3D()
            mol_aligned.copymol3D(mol_)
        theta0 += 5

    return mol

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

    if not len(anchoratoms) == 1:
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
            print('using DB distance of '+str(bondl))
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

## Loads ts M-L bond length from database and reports if compound is in DB
#  @param args Namespace of arguments
#  @param metal atom3D of atom 1 (usually a metal)
#  @param atom0 substrate connecting atom index
#  @param substrate Name of substrate
#  @param MLbonds M-L dictionary
#  @return Bond length in Angstroms, flag for exact DB match
def get_ts_MLdist_database(args,metal,subcatomsym,MLbonds):
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
    key.append((metal,oxs,spin,subcatomsym))
    key.append((metal,oxs,'-',subcatomsym))
    key.append((metal,'-',spin,subcatomsym))
    key.append((metal,'-','-',subcatomsym))
    found = False
    exact_match = False
    # search for data
    for kk in key:
        if (kk in MLbonds.keys()): # if exact key in dictionary
            bondl_sub = float(MLbonds[kk][0])
            bondl_m3D = float(MLbonds[kk][1])
            bondl_core3D = float(MLbonds[kk][2])
            found = True
            if (kk == ((metal,oxs,spin,subcatomsym))): ## exact match
               exact_match = True
            break
    if not found: # last resort covalent radii
        bondl_core3D = 1.8
        bondl_m3D = 1.25
        bondl_sub = 1.3
    if args.debug:
        print(MLbonds)
        print('key is ' + str(key))
    return bondl_core3D,bondl_m3D,bondl_sub,exact_match

## Loads ts M-L fsr from database and reports if compound is in DB
#  @param sym1 atomic symbol 1
#  @param sym2 atomic symbol 2
#  @param fsr_dict formal shortness ratio dictionary
#  @param ox oxidation state of the metal
#  @param spin spin multiplicity of the metal
#  @return fsr, flag for exact DB match
def get_ts_fsr_database(sym1, sym2, fsr_dict, ox=False, spin=False, bondl_m3D=False):
    # check for roman letters in oxstate
    if ox: # if defined put oxstate in keys
        if ox in romans.keys():
            oxs = romans[ox]
        else:
            oxs = ox
    else:
        oxs = '-'
    # check for spin multiplicity
    spin = spin if spin else '-'
    keys = []
    syms = [[sym1,sym2],[sym2,sym1]]
    oss = [[oxs,spin],[oxs,'-'],['-',spin],['-','-']]
    for sym in syms:
        for os in oss:
            key = sym + os
            keys.append(tuple(key))
    found = False
    exact_match = False
    # search for data
    for key in keys:
        if (key in fsr_dict.keys()): # if exact key in dictionary
            if bondl_m3D:
                fsr = float(fsr_dict[key][1])
            else:
                fsr = float(fsr_dict[key][0])
            found = True
            if (key == ((sym1, sym2, oxs, spin))): ## exact match
               exact_match = True
            break
    if not found: # last resort covalent radii
        fsr = 1
    return fsr, exact_match

# ## Loads ts M-L-S angle from database and reports if compound is in DB
# #  @param sym atomic symbol
# #  @param val remaining valency of the atom
# #  @param MLS_dict M-L-S angle dictionary
# #  @param ox oxidation state of the metal
# #  @param spin spin multiplicity of the metal
# #  @return bangle, flag for exact DB match
# def get_ts_MLSangle_database(sym,val,MLS_angle_dict,ox=False,spin=False):
#     # check for roman letters in oxstate
#     if ox: # if defined put oxstate in keys
#         if ox in romans.keys():
#             oxs = romans[ox]
#         else:
#             oxs = ox
#     else:
#         oxs = '-'
#     # check for spin multiplicity
#     spin = spin if spin else '-'
#     keys = []
#     syms = [[sym,str(val)]]
#     oss = [[oxs,spin],[oxs,'-'],['-',spin],['-','-']]
#     for sym_ in syms:
#         for os in oss:
#             key = sym_ + os
#             keys.append(tuple(key))
#     found = False
#     exact_match = False
#     # search for data
#     for key in keys:
#         if (key in MLS_angle_dict.keys()): # if exact key in dictionary
#             bangle = float(MLS_angle_dict[key][0])
#             found = True
#             if (key == ((sym,val,oxs,spin))): ## exact match
#                exact_match = True
#             break
#     if not found: # last resort covalent radii
#         bangle = 120
#     return bangle, exact_match

## Loads ts M-L-S angle from database and reports if compound is in DB
#  @param sym atomic symbol
#  @param sym valency
#  @param fsr_dict formal shortness ratio dictionary
#  @param ox oxidation state of the metal
#  @param spin spin multiplicity of the metal
#  @return fsr, flag for exact DB match
def get_ts_MLSangle_database(sym, val, MLS_dict, ox=False, spin=False):
    # check for roman letters in oxstate
    if ox: # if defined put oxstate in keys
        if ox in romans.keys():
            oxs = romans[ox]
        else:
            oxs = ox
    else:
        oxs = '-'
    # check for spin multiplicity
    spin = spin if spin else '-'
    keys = []
    syms = [[sym, val], [val, sym]]
    oss = [[oxs,spin],[oxs,'-'],['-',spin],['-','-']]
    for sym in syms:
        for os in oss:
            key = sym + os
            keys.append(tuple(key))
    found = False
    exact_match = False
    # search for data
    for key in keys:
        if (key in MLS_dict.keys()): # if exact key in dictionary
            bangle = float(MLS_dict[key][0])
            bdihedral = float(MLS_dict[key][1])
            found = True
            if (key == ((sym, val, oxs, spin))): ## exact match
               exact_match = True
            break
    if not found: # last resort covalent radii
        bangle = 120
        bdihedral = 90
    return bangle, bdihedral, exact_match

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

## Aligns a monodentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param substrate Name of substrate for dictionary lookup
#  @param sub3D mol3D of substrate
#  @param catoms List of substrate connecting atom indices
#  @param rempi Flag for pi-coordinating substrate
#  @param subpiatoms List of pi-coordinating atom indices in substrate
#  @param MLb Custom M-L bond length (if any)
##  @param ANN_flag Flag for ANN activation
##  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param i Ligand serial number
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand, updated list of M-L bond lengths
def align_sub(args,cpoints,core3D,coreref,sub3D,subcatoms,mligcatoms_ext,bangless,rempi,subpiatoms):
    bangle_m3Dsub = bangless[0][1]
# ,ANN_flag=False,ANN_bondl=[],this_diag=0,MLbonds=dict(),MLoptbds=[],i=0,EnableAutoLinearBend=True):
    corerefcoords = coreref.coords()
    # connection atom in sub3D
    atom0 = int(subcatoms[0])
    subcoords = sub3D.getAtom(atom0)
    if args.debug:
        print('atom0 is ' + str(atom0))
    # print substrate coordinates before translation
    if args.debug:
        print(corerefcoords)
        print(cpoints[0].coords())
        # print(atom0)
        # sub3D.printxyz()
    # translate ligand to overlap with backbone connecting point
    sub3D.alignmol(subcoords,cpoints[0])
    # determine bond length (database/cov rad/ANN)
    # bondl = get_MLdist(args,lig3D,atom0,ligand,coreref,MLb,i,ANN_flag,ANN_bondl,this_diag,MLbonds)
    # MLoptbds = []
    # bondl = 2
    # MLoptbds.append(bondl)
    # align ligand to correct M-L distance
    u = vecdiff(cpoints[0].coords(),corerefcoords)
    # sub3D = aligntoaxis2(sub3D, cpoint.coords(), corerefcoords, u, bondl)
    # print('cpoint is ' + str(cpoints[0].coords()) + '. corerefcoords is ' + str(corerefcoords) + '.')
    sub3D = aligntoaxis(sub3D, cpoints[0].coords(), corerefcoords, u)
    if args.debug:
        print('length of subpiatoms is ' + str(len(subpiatoms)))
    if rempi and len(subcatoms) > 1:
        # align linear (non-arom.) pi-coordinating ligand
        sub3D = align_linear_pi_sub(core3D,mligcatoms_ext,sub3D,atom0,subcatoms,bangle_m3Dsub)
        if args.debug:
            print('aligning a linear pi ligand at a bangle_m3Dsub of ' + str(bangle_m3Dsub))
    elif sub3D.natoms > 1:
        # align ligand center of symmetry
        sub3D = align_sub_firstbond(args,corerefcoords,sub3D,atom0,core3D,bangle_m3Dsub)
        print('bangle_m3Dsub is ' + str(bangle_m3Dsub))
        if sub3D.natoms > 2:
            # check for linear molecule and align
            sub3D = check_rotate_linear_lig(corerefcoords,sub3D,atom0)
            # check for symmetric molecule
            sub3D = check_rotate_symm_lig(corerefcoords,sub3D,atom0,core3D)
        # rotate around M-L axis to minimize steric repulsion
        # sub3D = rotate_MLaxis_minimize_steric_ts(mligcatoms_ext,sub3D,atom0,core3D)
        # rotate around L-Sub axis to minimize steric repulsion
        sub3D = rotate_MLaxis_minimize_steric(corerefcoords,sub3D,atom0,core3D)
    # rotate the substrate around the M-mlig axis to reduce repulsion
    mcoords = [core3D.getAtom(idx).coords() for idx in core3D.getBondedAtoms(mligcatoms_ext) if idx in core3D.findMetal()][0]
    rmref = vecdiff(mcoords, corerefcoords)
    sub3D_aligned = rotate_MLaxis_minimize_steric_ts(sub3D, subcoords, core3D, corerefcoords, rmref)
    # if args.debug:
    #     sub3D_aligned.printxyz()

    return sub3D_aligned

## Aligns an intramolecular monodentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param core3D mol3D of partially built complex
#  @param coreref atom3D of core reference atom
#  @param substrate Name of substrate for dictionary lookup
#  @param sub3D mol3D of substrate
#  @param catoms List of substrate connecting atom indices
#  @param rempi Flag for pi-coordinating substrate
#  @param subpiatoms List of pi-coordinating atom indices in substrate
#  @param MLb Custom M-L bond length (if any)
##  @param ANN_flag Flag for ANN activation
##  @param ANN_bondl ANN-predicted M-L bond length
#  @param this_diag ANN diagnostic object
#  @param MLbonds M-L bond dictionary
#  @param MLoptbds List of final M-L bond lengths
#  @param i Ligand serial number
#  @param EnableAutoLinearBend Flag for enabling automatic bending of linear ligands (e.g. superoxo)
#  @return mol3D of aligned ligand, updated list of M-L bond lengths
def align_intra_sub(args,core3D,subcatoms_ext,mligcatoms_ext,bondl_core3D,bondl_m3D,bondl_sub,bangless):
    bangle_m3D = bangless[0][0]
    bangle_m3Dsub = bangless[0][1]
    reacting_at_list = []
    reacting_at_list.append(subcatoms_ext)
    reacting_at_list.append(mligcatoms_ext)
    midx_list = core3D.findMetal()
    constr_list = []
    for midx in midx_list:
        constr_list.append(midx)
        for fidx in core3D.getBondedAtoms(midx):
            if fidx not in reacting_at_list:
                constr_list.append(fidx)
    sidx = core3D.getBondedAtoms(subcatoms_ext)[0]
    for i in core3D.getBondedAtoms(subcatoms_ext):
        if i in midx_list:
            sidx = i
    # bondl_sub_i = distance(core3D.getAtom(subcatoms_ext[0]).coords(),core3D.getAtom(sidx).coords())
    # bondl_m3D_i = distance(core3D.getAtom(subcatoms_ext[0]).coords(),core3D.getAtom(mligcatoms_ext[0]).coords())
    # bondl_core3D_i = distance(core3D.getAtom(mligcatoms_ext[0]).coords(),core3D.getAtom(midx_list[0]).coords())
    # bangle_m3D_i,u = rotation_params(core3D.getAtom(midx_list[0]).coords(),core3D.getAtom(mligcatoms_ext[0]).coords(),core3D.getAtom(subcatoms_ext[0]).coords())
    core3D.convert2OBMol()
    OBMol = core3D.OBMol
    # for bond in openbabel.OBMolBondIter(OBMol):
    #     idx1 = bond.GetBeginAtomIdx()
    #     idx2 = bond.GetEndAtomIdx()
    #     # if (idx1-1 in constr_list and idx2-1 in reacting_at_list) or (idx2-1 in constr_list and idx1-1 in reacting_at_list):
    #     #     bond.SetBO(0)
    #     if (idx1-1 == 31) and (idx2-1 == 33):
    #         bond.SetBO(2)

    # core3D.OBMol = OBMol
    # core3D.convert2mol3D()
    # core3D.convert2OBMol()
    # OBMol = core3D.OBMol
    ff = openbabel.OBForceField.FindForceField('mmff94')
    # for a_m3D in numpy.arange(bangle_m3D_i,bangle_m3D,-10).tolist():
    #     for b_sub in numpy.arange(bondl_sub_i,bondl_sub,-0.1).tolist():
    #         for b_m3D in numpy.arange(bondl_m3D_i,bondl_m3D,-0.1).tolist():
    #             for b_core3D in numpy.arange(bondl_core3D_i,bondl_core3D,-0.1).tolist():
    #                 constr = openbabel.OBFFConstraints()
    #                 for idx in constr_list:
    #                     constr.AddAtomConstraint(idx+1)
    #                 for subidx in subcatoms_ext:
    #                     constr.AddDistanceConstraint(subidx+1,sidx+1,b_sub) # bondl_sub
    #                     for mligcatoms in mligcatoms_ext:
    #                         constr.AddDistanceConstraint(subidx+1,mligcatoms+1,b_m3D) # bondl_m3D
    #                     for midx in midx_list:
    #                         constr.AddDistanceConstraint(midx+1,mligcatoms+1,b_core3D) # bondl_core3D
    #                         constr.AddAngleConstraint(midx+1,mligcatoms_ext[0]+1,subidx+1,a_m3D) # bangle_m3D
    constr = openbabel.OBFFConstraints()
    for idx in constr_list:
        constr.AddAtomConstraint(idx+1)
    subidx = subcatoms_ext
    constr.AddDistanceConstraint(subidx+1,sidx+1,bondl_sub) # bondl_sub
    if args.debug:
        print('Distance constraint between %s and %s is %s' %(subidx, sidx, bondl_sub))
    mligcatoms = mligcatoms_ext
    constr.AddDistanceConstraint(subidx+1,mligcatoms+1,bondl_m3D) # bondl_m3D
    if args.debug:
        print('Distance constraint between %s and %s is %s' % (subidx, mligcatoms, bondl_m3D))
    atnos = []
    for midx in midx_list:
        constr.AddDistanceConstraint(midx+1,mligcatoms+1,bondl_core3D) # bondl_core3D
        if args.debug:
            print('Distance constraint between %s and %s is %s' % (midx, mligcatoms, bondl_core3D))
        constr.AddAngleConstraint(midx+1,mligcatoms_ext+1,subidx+1,bangle_m3D) # bangle_m3D
        if args.debug:
            print('Angle constraint among %s, %s, and %s is %s' % (midx, mligcatoms_ext, subidx, bangle_m3D))
        atno = OBMol.GetAtom(midx+1).GetAtomicNum()
        atnos.append(atno)
        OBMol.GetAtom(midx+1).SetAtomicNum(14)
        # constr.AddDistanceConstraint(1,22,2.202)
    # printing obatom type
    # for obatom in openbabel.OBMolAtomIter(OBMol):
    #     print(obatom.GetType())
    s = ff.Setup(OBMol,constr)
    ff.SteepestDescent(500)
    for i, midx in enumerate(midx_list):
        atno = atnos[i]
        OBMol.GetAtom(midx+1).SetAtomicNum(atno)
    ff.GetCoordinates(OBMol)
    # printing obatom type
    # for obatom in openbabel.OBMolAtomIter(OBMol):
    #     print(obatom.GetType())
    core3D.OBMol = OBMol
    core3D.convert2mol3D()

    return core3D

## Aligns a bidentate substrate to core connecting atom coordinates
#  @param args Namespace of arguments
#  @param cpoint atom3D containing backbone connecting point
#  @param m3D mol3D of backbone template
#  @param core3D mol3D of partially built complex
#  @param mligcatom
#  @param sub3D mol3D of substrate
#  @param subcatoms List of substrate connecting atom indices
#  @param bangle_m3Dsub
#  @return mol3D of aligned ligand, updated lists of frozen atoms and M-L bond lengths
def align_dent2_sub(args,cpoints,m3D,core3D,mligcatom,sub3D,subcatoms,bangless):
    bangle_m3Dsub = bangless[0][1]
    corerefcoords = m3D.getAtom(mligcatom).coords()
    r0 = corerefcoords
    # get cis conformer by rotating rotatable bonds
    #lig3D = find_rotate_rotatable_bond(lig3D,catoms)
    # connection atom
    atom0 = subcatoms[0]
    atom1 = subcatoms[1]
    subcoords = sub3D.getAtom(atom1).coords()
    # translate ligand to match first connecting atom to backbone connecting point
    sub3D.alignmol(sub3D.getAtom(atom0),cpoints[0])
    # rotate the substrate to align with the cpoints
    r1 = sub3D.getAtom(atom0).coords() # atom0 coord
    r2 = sub3D.getAtom(atom1).coords()  # atom1 coord
    rx2 = cpoints[1].coords() # second cpoint coordxxx
    theta, u = rotation_params(r2, r1, rx2)
    sub3D = rotate_around_axis(sub3D, r1, u, 180-theta)
    # m3D_ = mol3D()
    # # for atom in m3D.atoms:
    # #     if atom.sym != 'X':
    # #         m3D_.addAtom(atom)
    # m3D_.copymol3D(m3D)
    # m3D_.copymol3D(sub3D)
    # rotate the sub3D along the two coordinating atoms to maximize the overlap between sub3D-Xs and the mcomplex
    r1 = sub3D.getAtom(atom0).coords() # atom0 coord
    r2 = sub3D.getAtom(atom1).coords()  # atom1 coord
    xidxes_m3D = sorted(m3D.findAtomsbySymbol('X'))
    xidxes_sub3D = sorted(sub3D.findAtomsbySymbol('X'))
    refidx = m3D.getBondedAtoms(xidxes_m3D[0])[0]
    refcoords = m3D.getAtomCoords(refidx)
    xcoords = sub3D.getAtomCoords(xidxes_sub3D[0])
    theta, u = rotation_params(refcoords, r1, xcoords)
    u = vecdiff(r1, r2)
    # print('theta is ' + str(theta))
    if theta < 90:
        theta = 180 - theta
    sub3D = rotate_around_axis(sub3D, r1, u, theta)

    sub3D_aligned = sub3D
    # sub3D_aligned = rotate_MLaxis_minimize_steric_ts(sub3D, subcoords, m3D, corerefcoords, rmref)
    # r1 = sub3D_aligned.getAtom(atom0).coords()  # atom0 coord
    # r2 = sub3D_aligned.getAtom(atom1).coords()  # atom1 coord
    # r12 = vecdiff(r2, r1)
    # subcoords = sub3D_aligned.getAtom(atom1).coords()
    # sub3D_aligned = rotate_MLaxis_minimize_steric_ts(sub3D_aligned, subcoords, m3D, subcoords, r12)

    return sub3D_aligned

## Crude rotations to improve alignment of the 2nd connecting atom of a bidentate substrate
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
    # tmp3D.writexyz('new')
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
    ANN_flag,ANN_bondl,ANN_reason,ANN_attributes, catalysis_flag = init_ANN(args,ligands,occs,dents,batslist,tcats,licores)



    this_diag.set_ANN(ANN_flag,ANN_reason,ANN_attributes, catalysis_flag)

    # freeze core
    for i in range(0,core3D.natoms):
        frozenats.append(i)

    # loop over ligands and begin functionalization
    # loop over ligands
    totlig = 0  # total number of ligands added
    ligsused = 0 # total number of ligands used
    loopcount = 0 # this counts the site occupations (I think?)
    subcatoms_ext = []
    mligcatoms_ext = []
    if args.mligcatoms:
        for i in range(len(args.mligcatoms)):
            mligcatoms_ext.append(0)
    initatoms = core3D.natoms # initial number of atoms in core3D
    if args.tsgen and (args.core.lower() in args.mlig):
        mligcatoms_ext[args.mlig.index(args.core.lower())] = 0
    for i,ligand in enumerate(ligands):
        if args.debug:
                print('************')
                print('loading ligand '+str(ligand) + ', number  ' + str(i) + ' of ' + str(len(ligands)))
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
        # Skip = False
        for j in range(0,occs[i]):
            if args.debug:
                print('loading copy '+str(j) + ' of ligand ' + ligand + ' with dent ' + str(dents[i]))
                print('totlig is ' + str(totlig))
                print('ligused is ' + str(ligsused))
                print('target BL is ' + str(ANN_bondl[ligsused]))
                print('******')
            denticity = dents[i]

            if not(ligand=='x' or ligand =='X') and (totlig-1+denticity < coord):
                # add atoms to connected atoms list
                catoms = lig.cat # connection atoms
                initatoms = core3D.natoms # initial number of atoms in core3D
                # if args.tsgen and (ligand in args.substrate) and (Skip is True):
                #     for k in args.subcatoms:
                #         subcatoms_ext.append(int(k) + initatoms)
                # if args.tsgen and (ligand.lower() in args.mlig) and (Skip is False):
                #     for k in args.mligcatoms:
                #         mligcatoms_ext.append(int(k) + initatoms)
                #     Skip = True
                if args.tsgen and (ligand.lower() in [substrate_i.lower() for substrate_i in args.substrate]):
                    for k in args.subcatoms:
                        subcatoms_ext.append(int(k) + initatoms)
                if args.tsgen and (ligand.lower() in [mlig_i.lower() for mlig_i in args.mlig]):
                    mligcatoms_ext[args.mlig.index(ligand.lower())] = int(args.mligcatoms[args.mlig.index(ligand.lower())]) + initatoms
                if args.debug:
                    print(ligand.lower())
                    print(args.mlig)
                    print(args.core.lower())
                    print('mligcatoms_ext is ' + str(mligcatoms_ext))
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
                if args.debug:
                    print('backbone atoms: ' + str(batoms))
                if (denticity == 1):
                    lig3D,MLoptbds = align_dent1_lig(args,cpoint,core3D,coreref,ligand,lig3D,catoms,rempi,ligpiatoms,MLb,ANN_flag,ANN_bondl[ligsused],this_diag,MLbonds,MLoptbds,i)
                    if args.debug:
                        print('adding monodentate at distance: ' + str(ANN_bondl[totlig]) + '/'+str(MLb)+ '/'+' at catoms ' + str(catoms))
                        print('printing ligand information')
                        print(lig3D.printxyz())
                elif (denticity == 2):
                    lig3D,frozenats,MLoptbds = align_dent2_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl[ligsused],this_diag,MLbonds,MLoptbds,frozenats,i)
                elif (denticity == 3):
                    lig3D,frozenats,MLoptbds = align_dent3_lig(args,cpoint,batoms,m3D,core3D,coreref,ligand,lig3D,catoms,MLb,ANN_flag,ANN_bondl[ligsused],this_diag,MLbonds,MLoptbds,frozenats,i)
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
                    #! TODO BCM defition slightly modified. Keep an eye for unexpected structures
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

                    bondl = get_MLdist(args,lig3D,atom0,ligand,m3D.getAtom(0),MLb,i,ANN_flag,ANN_bondl[ligsused],this_diag,MLbonds)
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
                        if args.debug:
                            print('a is not ff.lower, so adding atom:  ' + str(latdix+core3D.natoms)+  ' to freeze')
                        frozenats.append(latdix+core3D.natoms)
                # combine molecules
                core3D = core3D.combine(lig3D)
                core3D.convert2OBMol()
                core3D.convert2mol3D()
                # remove dummy cm atom if requested
                if rempi:
                    core3D.deleteatom(core3D.natoms-1)
                print('number of atoms in lig3D is ' + str(lig3D.natoms))
                if lig3D.natoms < 3:
                    frozenats += range(core3D.natoms-2,core3D.natoms)
                    print(str(range(core3D.natoms-2,core3D.natoms)) + ' are frozen.')
                if args.calccharge:
                    core3D.charge += lig3D.charge
                    print('core3D charge is ' + str(core3D.charge))
                # perform FF optimization if requested
                if args.debug:
                    print('saving a copy of the complex named complex_'+str(i)+'_'+str(j) + '.xyz')
                    core3D.writexyz('complex_'+str(i)+'_'+str(j) + '.xyz')

                if 'a' in args.ffoption:
                    if args.debug:
                        print('FF optimizing molecule after placing ligand')
                        print('in the relax function, passing connected atoms list: ' + str(connected))
                    #(ff,mol,connected,constopt,frozenats,frozenangles,mlbonds,nsteps,debug=False):
                    core3D,enc = ffopt(ff=args.ff,\
                                        mol=core3D,\
                                        connected=connected,\
                                        constopt=1,\
                                        frozenats=frozenats,\
                                        frozenangles=freezeangles,\
                                        mlbonds=MLoptbds,\
                                        nsteps='Adaptive',\
                                        debug=args.debug)
                    if args.debug:
                        print('saving a copy of the complex named complex_'+str(i)+'_'+str(j) + '_ff.xyz')
                        core3D.writexyz('complex_'+str(i)+'_'+str(j) + '_ff.xyz')
                if args.debug:
                    print('done with pair of inds '+str(i)+' and '+str(j))
                    print('**************************')
            totlig += denticity
            ligsused += 1
    # perform FF optimization if requested
    if 'a' in args.ffoption:
        print('Performing final FF opt')
        # idxes
        midx = core3D.findMetal()[0]
        fidxes = core3D.getBondedAtoms(midx)
        fidxes = [fidx for fidx in fidxes if core3D.getAtom(fidx).sym != 'H']
        fsyms = [core3D.getAtom(fidx).sym for fidx in fidxes]
        print(fsyms)
        # constraints
        constr = openbabel.OBFFConstraints()
        for idx in fidxes + [midx]:
            constr.AddAtomConstraint(idx+1)
        # ff
        ff = openbabel.OBForceField.FindType('UFF')
        core3D.convert2OBMol()
        obmol = core3D.OBMol
        flag = ff.Setup(obmol,constr)
        if flag:
            ff.SteepestDescent(5000)
        ff.GetCoordinates(obmol)
        core3D.OBMol = obmol
        core3D.convert2mol3D()

        # core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
    return core3D,complex3D,emsg,this_diag,subcatoms_ext,mligcatoms_ext

## Main substrate placement routine for transition state
#  @param args Namespace of arguments
#  @oaram core3D core3D of the mcomplex
#  @param substrate a list of substrates
#  @param sub_i the index of the substrate
#  @param subcatoms the atom index in the substrate where the mcomplex is attached to
#  @param mlig the
#  @param subcatoms_ext the atom index, compounded with the number of atoms in the mcomplex,
#  in the substrate where the mcomplex is attahced to
#  @param mligcatoms_ext the atom index, compounded with the number of atoms in the mcomplex,
#  in the mcoomplex where the substrate is attahced to
#  @return mol3D of built complex, list of all mol3D ligands and core, error messages
def msubcomplex(args, core3D, substrate, sub_i, subcatoms, mlig, subcatoms_ext, mligcatoms_ext):

    globs = globalvars()
    subcores = getsubcores()
    this_diag = run_diag()
    if globs.debug:
        print '\nGenerating TS complex with substrate and mlig:',substrate,mlig
    if args.gui:
        args.gui.iWtxt.setText('\nGenerating complex with core:'+args.core+' and ligands: '+ ' '.join(arg.ligs)+'\n'+
                               args.gui.iWtxt.toPlainText())
        args.gui.app.processEvents()
    # import gui options
    if args.gui:
        from Classes.mWidgets import mQDialogWarn
    # initialize variables
    emsg, complex3D = False, [] # function returns
    occs = 1        # currently only one substrate with the occurance of one is supported.
    catsmi = []     # SMILES substrates connection atoms
    smilessub = 0   # count how many smiles strings
    cats0 = []      # connection atoms for each substrate
    dentl = []      # denticity of substrates
    # tcats = []      # list of connection atoms for all substrates
    connected = []  # indices in core3D of substrate atoms connected to metal
    frozenats, freezeangles = [], False  # ff frozen variables
    MLoptbds = []   # list of bond lengths
    rempi = False   # remove dummy pi orbital center of mass atom
    backbatoms, batslist, bats = [], [], [] # backbond atoms variables
    bondls_m3D, bangles_m3D = [], [] # molecular parameters
    charge = 0
    if args.calccharge:
        charge = core3D.charge
    elif args.charge:
        charge = args.charge[0]
    ox = args.oxstate if args.oxstate else False
    spin = args.spin if args.spin else False
    # load substrate
    sub, subcatoms, emsg = substr_load(substrate[0], sub_i, subcatoms)
    sub.convert2mol3D()
    # # remove dummy atoms in core3D
    for atom_i, atom in enumerate(core3D.atoms):
        asym = atom.sym
        if asym == 'X':
            core3D.deleteatom(atom_i)
    # obtain rxn types
    if len(sub.grps) > 0:
        rxn_types = sub.grps
    else:
        rxn_types = ['inter'] * len(subcatoms)
    # obtaining molecular parameters by rxn types
    bondlss = []
    bangless = []
    bdihedralss = []
    if len(rxn_types) != len(mligcatoms_ext):
        mligcatoms_ext = mligcatoms_ext * len(rxn_types)
    for rxn_type_i, rxn_type in enumerate(rxn_types):
        #
        ## load bond data
        #
        # load dicts
        fsr_dict = loaddata_ts('/Data/ML_FSR_for_' + rxn_type + '.dat')
        MLS_dict = loaddata_ts('/Data/MLS_FSR_for_' + rxn_type + '.dat')
        # get indexes
        midxes = core3D.findMetal()  # list of metal indexes
        mligcatom_ext = mligcatoms_ext[rxn_type_i]
        mcoord = core3D.getAtom(midxes[0]).coords()
        d0 = 10
        for i, fidx in enumerate(core3D.getBondedAtoms(mligcatom_ext)):
            fcoord = core3D.getAtom(fidx).coords()
            d = distance(fcoord, fcoord)
            if d < d0:
                d0 = d
                mligfidx_ext = fidx
        subcatom = int(subcatoms[rxn_type_i])
        bnum0 = 4
        d0 = 10
        subfidx = 0
        for i, fidx in enumerate(sub.getBondedAtoms(subcatom)):
            if fidx in subcatoms:
                subfidx = fidx
                break
            else:
                if len(subcatoms) > rxn_type_i + 1:
                    d = distance(sub.getAtomCoords(fidx), sub.getAtomCoords(subcatoms[rxn_type_i+1]))
                    if d < d0:
                        d0 = d
                        subfidx = fidx
                else:
                    bnum = len(sub.getBondedAtoms(fidx))
                    if bnum < bnum0:
                        bnum0 = bnum
                        subfidx = fidx
        #! TODO this part should be contracted by using the same method
        # get atom symbols
        msym = core3D.getAtom(midxes[0]).sym
        mligcatomsym = core3D.getAtom(mligcatom_ext).sym
        mligfidxsym = core3D.getAtom(mligfidx_ext).sym
        subcatomsym = sub.getAtom(subcatom).sym
        subfidxsym = sub.getAtom(subfidx).sym
        # get valence electron counts
        mligcatomve = int(amassdict[mligcatomsym][3])
        mligfidxve = int(amassdict[mligfidxsym][3])
        subcatomve = int(amassdict[subcatomsym][3])
        subfidxve = int(amassdict[subfidxsym][3])
        # get the number of bonded atoms
        mligcatombnum = len(core3D.getBondedAtoms(mligcatom_ext))
        subcatombnum = len(sub.getBondedAtoms(subcatom))
        mligfidxbnum = len(core3D.getBondedAtoms(mligfidx_ext))
        subfidxbnum = len(sub.getBondedAtoms(subfidx))
        # get remaining valency
        octet = 8 if mligcatomsym != 'H' else 2
        mligcatomval = (octet - mligcatomve - mligcatombnum)
        octet = 8 if subcatomsym != 'H' else 2
        subcatomval = (octet - subcatomve - subcatombnum)
        if mligfidxsym in ['V', 'Cr', 'Mn', 'Fe', 'Co']:
            octet = 18
            multiplier = 2
        elif mligfidxsym in ['Ni', 'Cu']:
            octet = 14
            multiplier = 2
        elif mligfidxsym != 'H':
            octet = 8
            multiplier = 1
        else:
            octet = 2
            multiplier = 1
        mligfidxval = (octet - mligfidxve - mligfidxbnum * multiplier)
        octet = 8 if subfidxsym != 'H' else 2
        subfidxval = (octet - subfidxve - subfidxbnum)
        ## preping for bondls
        # get covalent radius
        mr = float(amassdict[msym][2])
        mligcatomr = float(amassdict[mligcatomsym][2])
        mligfidxr = float(amassdict[mligfidxsym][2])
        subcatomr = float(amassdict[subcatomsym][2])
        subfidxr = float(amassdict[subfidxsym][2])
        # get sum of cov radii
        sumr_m_mlig = mr + mligcatomr
        sumr_mlig_mligfidx = mligcatomr + mligfidxr
        sumr_mlig_sub = mligcatomr + subcatomr
        sumr_sub_subfidx = subcatomr + subfidxr
        # keys for bondls, in the order of bondl_core3D, bondl_m3D, and bondl_sub3D
        keys = [[mligcatomsym + str(mligcatomval), mligfidxsym + str(mligfidxval), sumr_mlig_mligfidx],
                [mligcatomsym + str(mligcatomval), subcatomsym + str(subcatomval), sumr_mlig_sub],
                [subcatomsym + str(subcatomval), subfidxsym + str(subfidxval), sumr_sub_subfidx]]
        # obtain bondls and bangles
        # default ratio
        fsr = 1
        bondls = []
        for i, key in enumerate(keys):
            bondl_m3D = False
            sym1 = key[0]
            sym2 = key[1]
            sumr = key[2]
            if i == 1:
                bondl_m3D = True
            fsr, exact_match = get_ts_fsr_database(sym1, sym2, fsr_dict, ox, spin, bondl_m3D=bondl_m3D)
            bondl = fsr * sumr
            bondls.append(bondl)
        bondlss.append(bondls)
        print('keys are ' + str(keys))
        print('bondlss are ' + str(bondlss))
        ## preping for bangles
        # keys for bangles, in  the order of bangle_m3D, bangle_m3Dsub, bangle_core3D, and bangle_sub
        # keys = [[mligcatomsym, mligcatomval], [subcatomsym, subcatomval], [mligfidxsym, mligfidxval],
        #         [subfidxsym, subfidxval]]
        keys = [[mligcatomsym, str(mligcatomval)], [subcatomsym, str(subcatomval)], [mligfidxsym, str(mligfidxval)],
                [subfidxsym, str(subfidxval)]]
        bdihedrals = []
        bangles = []
        for key in keys:
            sym = key[0]
            val = key[1]
            bangle, bdihedral, exact_match = get_ts_MLSangle_database(sym, val, MLS_dict, ox, spin)
            bangles.append(bangle)
            bdihedrals.append(bdihedral)
        bangless.append(bangles)
        bdihedralss.append(bdihedrals)
        # bangless = [[120, 180, 180, 109]]
        print('keys are ' + str(keys))
        print('bangless are ' + str(bangless))
        print('bdihedralss are ' + str(bdihedralss))
    # determine if the sub in lig
    sub_in_lig = False
    if args.lig:
        if (args.substrate[0].lower() in [lig_i.lower() for lig_i in args.lig]):
            sub_in_lig = True
    if not sub_in_lig:
        # freeze core
        for i in range(0,core3D.natoms):
            frozenats.append(i)
        # freeze key atoms in substrate
        if args.debug:
            print('subcatoms after init_substrate is ' + str(subcatoms))
        if len(subcatoms) > 1:
            for subcatom in subcatoms:
                if isinstance(subcatom,int):
                    frozenats.append(core3D.natoms + subcatom)
        else:
            frozenats.append(core3D.natoms + int(subcatoms[0]))
            for bondedat in sub.getBondedAtoms(int(subcatoms[0])):
                frozenats.append(core3D.natoms + bondedat)
        # compute number of connecting points required
        cpoints_required = len(bondlss)
        # load core and initialize template
        # also adjust the specified distance in the core3D and distance and angle in the connection point
        m3D, core3D, geom, backbatoms, coord, corerefatoms = init_mcomplex_template(args, core3D, cpoints_required,
                                                                                    mligcatoms_ext, bondlss,bangless,
                                                                                    bdihedralss)

        totsub = 0  # total number of substrates added
        subsused = 0
        ## initialize ligand
        sub3D, rempi, subcatoms, subpiatoms = init_substrate(args, sub, subcatoms, bondlss, bangless, bdihedralss)
        for subcatom in subcatoms:
            frozenats.append(core3D.natoms + subcatom)
            for batidx in sub3D.getBondedAtoms(subcatom):
                frozenats.append(core3D.natoms + batidx)
        if emsg:
            return False, emsg
        for j in range(0, occs):
            denticity = sub.denticity
            # print('denticity is ' + str(denticity))
            if not(substrate == 'x' or substrate == 'X'):
                # add atoms to connected atoms list
                # catoms = sub.cat # connection atoms
                # initatoms = core3D.natoms # initial number of atoms in core3D
                # for at in catoms:
                #     connected.append(initatoms+at)
                # initialize variables
                mcoords = core3D.getAtom(midxes[0]).coords() # metal c oordinates in backbone
                atom0, r0, r1, r2, r3 = 0, mcoords, 0, 0, 0 # initialize variables
                coreref = corerefatoms.getAtom(totsub)
                # connecting point in backbone to align ligand to
                # batoms = get_batoms(args,batslist,subsused)
                cpoints = m3D.atoms[m3D.natoms-denticity:]
                # attach ligand depending on the denticity
                # optimize geometry by minimizing steric effects
                if (denticity == 1):
                    sub3D = align_sub(args, cpoints, core3D, coreref, sub3D, subcatoms, mligcatoms_ext[0], bangless,
                                      rempi, subpiatoms)
                elif (denticity == 2):
                    # for mligcatom_ext in mligcatoms_ext:
                    #     if mligcatom_ext not in core3D.findMetal():
                    #         break
                    # print('cpoints are ' + str([cpoint.coords() for cpoint in cpoints]))
                    # print('mligcatom_ext is ' + str(mligcatom_ext))
                    # print('subcatoms are ' + str(subcatoms))
                    sub3D = align_dent2_sub(args, cpoints, m3D, core3D, mligcatom_ext, sub3D, subcatoms, bangless)
                auxm = mol3D()
                auxm.copymol3D(sub3D)
                complex3D.append(auxm)
                # if 'a' not in sub.ffopt.lower():
                #     for latdix in range(0,sub3D.natoms):
                #         frozenats.append(latdix+core3D.natoms)
                # core3D_ = mol3D()
                # core3D_.copymol3D(core3D)
                # core3D_ = core3D_.combine(sub3D)
                # core3D_.convert2mol3D()
                # remove dummy cm atom if requested
                # if rempi:
                #     core3D.deleteatom(core3D.natoms-1)
                ## final bondl adjustments
                # loading model
                # alpha = 0.005408009919423026
                # gamma = 0.017841136136665325
                # stat_dict, impt_dict, train_dict, test_dict, perm_dict, regr = gbr_model_training(csvf, 1, 2, 2)
                # stat_dict, impt_dict, train_dict, test_dict, perm_dict, regr = krr_model_training(csvf, 1, 2, alpha=alpha, gamma=gamma)
                # bondl prediction
                # spin = int(args.spin)
                # # bondl_dict, ds = ML_model_predict(core3D, spin, train_dict, stat_dict, impt_dict, regr)
                # bondl_dict, bondl_m3D, ds1, ds2 = krr_model_predict(core3D_, spin, mligcatoms_ext[0])
                if args.cdxml:
                    # combine molecules into a crude msubcomplex
                    core3D_ = mol3D()
                    core3D_.copymol3D(core3D)
                    liglist, ligdents, ligcons = ligand_breakdown(core3D_)
                    ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
                        ax_con_list, eq_con_list, built_ligand_list = ligand_assign(core3D_, liglist, ligdents, ligcons)
                    if globs.testTF():
                        ## new RACs-ANN
                        from molSimplify.Scripts.tf_nn_prep import invoke_ANNs_from_mol3d
                    result_dict = invoke_ANNs_from_mol3d(core3D_, oxidation_state=3) # result_dict, {"ls_bl": r_ls, "hs_bl": r_hs, "split": split, "distance": split_dist}
                    # for ds_i, ds in enumerate(ds1):
                    #     print('The euclidean distance between the built M-L' + str(ds_i) + ' and the training data is ' + str(min(ds)))
                    eu_dist = result_dict['distance']
                    print('The euclidean distance between the built rOH and the training data is ' + str(eu_dist))
                    midx = core3D_.findMetal()[0]
                    bondls = result_dict['ls_bl'][0]
                    for bondl_i, bondl in enumerate(bondls):
                        if bondl_i < 2:
                            idxes = ax_con_list[bondl_i]
                            for idx in idxes:
                                core3D_.BCM(idx, midx, bondl)
                                # print('Adjusting the distance between %s and %s to %s' %(idx, midx, bondl))
                        else:
                            idxes = eq_con_list[bondl_i-2]
                            for idx in idxes:
                                core3D_.BCM(idx, midx, bondl)
                                # print('Adjusting the distance between %s and %s to %s' % (idx, midx, bondl))
                    # # refined msubcomplex build
                    m3D, core3D, geom, backbatoms, coord, corerefatoms = init_mcomplex_template(args, core3D_,
                                                                                                cpoints_required,
                                                                                                mligcatoms_ext,
                                                                                                bondlss,
                                                                                                bangless)
                    # m3D.writexyz(fpath + '/init_mcomplex.xyz')

                    if (denticity == 1):
                        sub3D = align_sub(args, cpoints, core3D, coreref, sub3D, subcatoms, mligcatoms_ext[0], bangless,
                                          rempi, subpiatoms)

                    elif (denticity == 2):
                        # for mligcatom_ext in mligcatoms_ext:
                        #     if mligcatom_ext not in core3D_.findMetal():
                        #         break
                        # print('mligcatom_ext is ' + str(mligcatom_ext))
                        sub3D = align_dent2_sub(args, cpoints, m3D, core3D, mligcatom_ext, sub3D, subcatoms, bangless)
                auxm = mol3D()
                auxm.copymol3D(sub3D)
                complex3D.append(auxm)
                # if 'a' not in sub.ffopt.lower():
                #     for latdix in range(0,sub3D.natoms):
                #         frozenats.append(latdix+core3D.natoms)
                # combine molecules into a crude msubcomplex
                core3D = core3D.combine(sub3D)
                core3D.convert2mol3D()
                # core3D.writexyz(fpath + '/msubcomplex.xyz')
                # remove dummy cm atom if requested
                # if rempi:
                #     core3D.deleteatom(core3D.natoms-1)
                if args.calccharge:
                    core3D.charge += sub3D.charge
                totsub += denticity
                subsused += 1
        # perform FF optimization if requested
        if 'a' in args.ffoption:
            print('Performing final FF constrained opt for the metal-substrate complex')
            connected = []
            freezeangles = []
            core3D, enc = ffopt(args.ff, core3D, connected, 1, frozenats, freezeangles, MLoptbds, 'Adaptive',
                                args.debug)
            # core3D.writexyz(fpath + '/msubcomplex_ff.xyz')
    else:
        core3D = align_intra_sub(args, core3D, subcatoms_ext[0], mligcatoms_ext[0], bondlss, bangless)
        ## reserved for future conformer search development
        # for cycle in range(5):
        #     if 'a' in args.ffoption:
        #         print('Performing five conformer search: ' + str(cycle) + ' out of 5.')
        #         connected = []
        #         freezeangles = []
        #         core3D,enc = conformer_search(args.ff,core3D,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
        #         # core3D,enc = ffopt(args.ff,core3D,connected,1,frozenats,freezeangles,MLoptbds,'Adaptive',args.debug)
    return core3D, complex3D, subcatoms, emsg, this_diag

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
    subcores = getsubcores()
    # build structure
    sanity = False
    this_diag = run_diag()
    if (ligands):
        core3D,complex3D,emsg,this_diag,subcatoms_ext,mligcatoms_ext = mcomplex(args,ligands,ligoc,licores,globs)
        if args.debug:
            print('subcatoms_ext are ' + str(subcatoms_ext))
        name_core = args.core
        if emsg:
            return False,emsg
    else:
        print('You specified no ligands. The whole mcomplex is read from the core.')
        # read mol3D from core
        core3D = mol3D()
        if '.xyz' in args.core:
            core3D.readfromxyz(args.core)
        else:
            atom = atom3D()
            atom.__init__(Sym=args.core, xyz=[0, 0, 0])
            core3D.addAtom(atom)
        name_core = args.core
        if args.tsgen:
            if args.subcatoms:
                subcatoms_ext = core3D.natoms + int(args.subcatoms[0])
            else:
                args.subcatoms = [0]
                subcatoms_ext = core3D.natoms + 0
            midxes = core3D.findMetal()
            mligcatoms_ext = []
            fidxes = core3D.getBondedAtoms(midxes[0])
            for fidx in fidxes:
                atom = core3D.getAtom(fidx)
                asym = atom.sym
                if asym == 'O':
                    sidxes = core3D.getBondedAtoms(fidx)
                    if len(sidxes) == 1:
                        mligcatoms_ext.append(fidx)
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
    fname = name_complex(rootdir,name_core,args.geometry,ligands,ligoc,sernum,args,nconf,sanity)
    print('fname is ' + fname)
    # generate ts 
    if (args.tsgen):
        substrate = []
        for i in args.substrate:
            substrate.append(i)
        subcatoms = False
        subcatoms_all = ['0']
        if args.subcatoms:
            if 'hat' in [str(subcatom).lower() for subcatom in args.subcatoms]:
                sub,emsg,subcatoms = substr_load(substrate[0],0,subcatoms)
                sub.convert2mol3D()
                subcatoms_all = sub.getHs()
                if args.debug:
                    print(subcatoms_all)
            if 'epo' in [str(subcatom).lower() for subcatom in args.subcatoms]:
                sub,emsg,subcatoms = substr_load(substrate[0],0,subcatoms)
                sub.convert2mol3D()
                subcatoms_all = sub.getC2Cs()
                if args.debug:
                    print(subcatoms_all)
            else:
                subcatoms = []
                for i in args.subcatoms:
                    subcatoms.append(int(i))
        mlig = args.mlig
        mligcatoms = args.mligcatoms
        if not args.subcatoms and len([i for i in subcores.keys() if i.subname == substrate[0]]) > 0:
            nruns = len([i for i in subcores.keys() if i.subname == substrate[0]]) > 0
        else:
            nruns = len(subcatoms_all)
        for runs in range(nruns):
            if 'hat' in args.subcatoms or 'epo' in args.subcatoms:
                sub_i = 0
                if type(subcatoms_all[runs]) == list:
                    subcatoms = subcatoms_all[runs]
                else:
                    subcatoms = [subcatoms_all[runs]]
                if args.debug:
                    print('the subcatoms for run ' + str(runs) + ' is ' + str(subcatoms))
            else:
                sub_i = runs
            # if args.conformer:
            #     ncycles = 5
            # else:
            #     ncycles = 1
            # for cycles in range(ncycles):
            core3D_i = mol3D()
            core3D_i.copymol3D(core3D)
            core3D_i,complex3D,subcatoms,emsg,this_diag = msubcomplex(args,core3D_i,substrate,sub_i,subcatoms,mlig,subcatoms_ext,mligcatoms_ext)
            fname = name_ts_complex(rootdir,name_core,args.geometry,ligands,ligoc,substrate,subcatoms,mlig,mligcatoms,sernum,args,nconf,sanity)
            if args.debug:
                print('fname is ' + str(fname))
            # write xyz file
            core3D_i.writexyz(fname)
            strfiles.append(fname)
            # write report file
            this_diag.set_mol(core3D_i)
            this_diag.write_report(fname+'.report')
            # write input file from command line arguments
            getinputargs(args,fname)
        return strfiles, emsg, this_diag
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

    #for ligs in ligands:
        #print(ligs)
        #sad

    if args.smicat: #or True:
        if sum([len(i)>1 for i in args.smicat]) > 0:
        #if True:
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
    elif args.tsgen:
        Nogeom = len(args.subcatoms)
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
