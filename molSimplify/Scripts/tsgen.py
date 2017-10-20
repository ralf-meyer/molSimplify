#!/usr/bin/env python
# Written by Terry Gani for HJK Group
# Dpt of Chemical Engineering, MIT

# This script generates TS guesses based on a user-specified complex (core), substrate and reaction.

# Currently, only the hydrogen abstraction reaction is supported. 

# import custom modules
from molSimplify.Scripts.geometry import *
from molSimplify.Scripts.io import *
from molSimplify.Scripts.nn_prep import *
from molSimplify.Classes.globalvars import *
from molSimplify.Classes.rundiag import *
from molSimplify.Classes import globalvars
from molSimplify.Classes import mol3D
from molSimplify.Informatics.decoration_manager import*
# import standard modules
import os, sys
from pybel import *
import openbabel
from pkg_resources import resource_filename, Requirement
import random, itertools
from numpy import log, arccos, cross, dot, pi

angopt = 180

def getconnection2(core,catom,Midx,BL,angopt):
    angcoeff = 0.1 # coefficient of MXY term in connecting point optimization
    ### add fake atoms for catoms
    ncore = core.natoms
    # add fake atom in local centermass axis
    Ocoords = core.getAtom(catom).coords()
    Mcoords = core.getAtom(Midx).coords()
    backbcoords = alignPtoaxis(Ocoords,Ocoords,vecdiff(Ocoords,Mcoords),BL)
    # manually find best positioning
    am = mol3D()
    am.addAtom(atom3D('C',backbcoords))
    corerem = mol3D()
    corerem.copymol3D(core)
    corerem.deleteatom(catom) # complex less O atom
    setopt = []
    mdist = -1
    for itheta in range(1,359,5):
        for iphi in range(1,179,5):
            P = PointTranslateSph(Ocoords,backbcoords,[BL,itheta,iphi])
            am.getAtom(0).setcoords(P)
            ang = 180-vecangle(vecdiff(Ocoords,Mcoords),vecdiff(P,Ocoords))
            if corerem.mindist(am) < 1: # reject if too close
                d0 = -2
            else:
                d0 = distance(core.centermass(),am.getAtomCoords(0))+0.5*log(corerem.mindist(am)-1)-angcoeff*(ang-angopt)*(ang-angopt)
            if d0 > mdist:
                mdist = d0
                setopt = am.coordsvect()
    core.addAtom(atom3D('C',setopt[0])) # Position placeholder, atom is deleted later
    connPts = core.getAtom(ncore).coords()
    return connPts

def getconnection1(core,Midx,BL):
	# get connecting point for oxidative addition based on a trained/fixed bond length and maximizing the smallest angle with other connecting atoms
	corecm = core.centermass()
	ncore = core.natoms
	nligs = core.getBondedAtoms(Midx)
	Mcoords = core.getAtom(Midx).coords()
	# manually find best positioning
	cpoint = []
	mdist = -1
	for itheta in range(1,359,5):
		for iphi in range(1,179,5):
			P = PointTranslateSph(Mcoords,Mcoords,[BL,itheta,iphi])
			d = []
			for iligs in nligs:
				d.append(vecangle(vecdiff(Mcoords,P),vecdiff(Mcoords,core.getAtomCoords(iligs))))
			d0 = min(d)
			if d0 > mdist:
				mdist = d0
				cpoint = P
	return cpoint

def distort_substr(substr,substreact,XYBL):
    adjsidx = substr.getBondedAtoms(substreact)[0]
    substr.BCM(substreact,adjsidx,XYBL)
    return substr    

def tsgen(mode,args,rootdir,core,substr,compreact,substreact,globs):
    # INPUT
    #   - mode: TS generation mode (see rungen.py)
    #   - args: placeholder for input arguments
    #   - rootdir: directory of current run
    #   - core: core mol3D (generated in rungen.py)
    #   - substr: substrate mol3D
    #   - compreact: complex reacting atom(s) (see rungen)
    #   - substreact: substrate reacting atom(s) (see rungen)
    #   - globs: class with global variables
    # OUTPUT
    #   - strfiles: list of xyz files generated
    #   - emsg: error messages
    emsg = False
    this_diag = run_diag()
    ############ LOAD DICTIONARIES ############
    mcores = getmcores()
    licores = getlicores()
    bindcores = getbcores()
    ########## END LOAD DICTIONARIES ##########
    strfiles = []
    ########## START FUNCTIONALIZING ##########

    sanity = False
    if mode == 2:
        emsg = 'Sorry, this mode is not supported yet. Exiting...'
        return strfiles, emsg, this_diag
    elif mode == 1: # oxidative addition of a single group
		# get first connecting point
        MLBL = 0.9*(core.getAtom(compreact).rad + substr.getAtom(substreact).rad)
        # In future this will be trained. Using sum cov rad for now.
        cpoint = getconnection1(core,compreact,MLBL)
        # distort substrate molecule
        adjsidx = substr.getBondedAtoms(substreact)[0]
        XYBL = 1.1*(substr.getAtom(substreact).rad + substr.getAtom(adjsidx).rad)
        substr = distort_substr(substr,substreact,XYBL)
        # align substrate molecule
        MXYang = 135
        substr.alignmol(substr.getAtom(substreact),atom3D('H',cpoint))
        tmp3D = mol3D()
        tmp3D.copymol3D(core)
        tmp3D.addAtom(atom3D('Cl',cpoint))
        ligalignpt = getconnection2(tmp3D,tmp3D.natoms-1,compreact,XYBL,MXYang)
        theta,u = rotation_params(substr.getAtomCoords(adjsidx),cpoint,ligalignpt)
        substr = rotate_around_axis(substr,cpoint,u,180-theta)
        # for substrates with more than 2 atoms, rotate around X-Y axis to minimize steric repulsion
        if substr.natoms > 2:
            dtheta = 2
            optmax = -9999
            totiters = 0
            substrb = mol3D()
            substrb.copymol3D(substr)
            # check for minimum distance between atoms and center of mass distance
            while totiters < 180:
                substr = rotate_around_axis(substr,cpoint,vecdiff(cpoint,ligalignpt),dtheta)
                d0 = substr.mindist(tmp3D) # try to maximize minimum atoms distance
                d0cm = substr.distance(tmp3D) # try to maximize center of mass distance
                iteropt = d0cm+10*log(d0) # optimization function
                if (iteropt > optmax): # if better conformation, keep
                    substrb = mol3D()
                    substrb.copymol3D(substr)
                optmax = iteropt
                totiters += 1
            substr = substrb
    elif mode == 3: # abstraction
        # distort A-B bond
        adjcidx = core.getBondedAtoms(compreact)[0]
        ABBL = 1.1*(core.getAtom(compreact).rad + core.getAtom(adjcidx).rad)
        core.BCM(compreact,adjcidx,ABBL)
        # set B-X distance 
        BXBL = 1.1*(substr.getAtom(substreact).rad + core.getAtom(compreact).rad)
        Bcoords = core.getAtom(compreact).coords()
        # get connecting point of abstracted atom
        angopt = 180
        cpoint = getconnection2(core,compreact,adjcidx,BXBL,angopt)
        core.deleteatom(core.natoms-1) # delete added atom
        # copy substrate
        substr_copy = mol3D()
        substr_copy.copymol3D(substr)
        # distort substrate molecule
        adjsidx = substr.getBondedAtoms(substreact)[0]
        XYBL = 1.1*(substr.getAtom(substreact).rad + substr.getAtom(adjsidx).rad)
        substr = distort_substr(substr,substreact,XYBL)
        # align substrate according to connection atom and shadow atom
        substr.alignmol(substr.getAtom(substreact),atom3D('H',cpoint))
        # perform rotations
        mcoords = core.getAtomCoords(compreact)
        r0 = core.getAtomCoords(compreact)
        r1 = substr.getAtomCoords(substreact)
        r2 = substr.centermass() # center of mass
        if not r2:
            emsg = 'Center of mass calculation for substrate failed. Check input.'
            print emsg
        at = substr_copy.getBondedAtoms(substreact)
        rrot = r1
        theta,u = rotation_params(r0,r1,r2)
        # align center of mass of local environment
        if (substr.natoms > 1):
            substrb = mol3D()
            substrb.copymol3D(substr)
            ####################################
            # center of mass of local environment (to avoid bad placement of bulky ligands)
            auxmol = mol3D()
            for i in at:
                auxmol.addAtom(substr.getAtom(i))   
            r2 = auxmol.centermass() # overwrite global with local centermass
            theta,u = rotation_params(r0,r1,r2)
            ####################################
            # rotate around axis and get both images
            substr = rotate_around_axis(substr,rrot,u,theta)
            substrb = rotate_around_axis(substrb,rrot,u,theta-180)
            d2 = distance(mcoords,substr.centermass())
            d1 = distance(mcoords,substrb.centermass())
            substr = substr if (d1 < d2) else substrb # pick better one
        if substr.natoms > 2:
            #####################################
            # check for linear molecule
            auxm = mol3D()
            for i in at: # BUGGY, SHOULD BE CHECKED IN THE PRE-DISTORTED SUBSTRATE INSTEAD
                auxm.addAtom(substr.getAtom(i))
            if auxm.natoms > 1:
                r0 = substr.getAtom(substreact).coords()
                r1 = auxm.getAtom(0).coords()
                r2 = auxm.getAtom(1).coords()
                if checkcolinear(r1,r0,r2):
                    # we will rotate so that angle is right
                    theta,urot = rotation_params(r1,mcoords,r2)
                    theta = vecangle(vecdiff(r0,mcoords),urot)
                    substr = rotate_around_axis(substr,r0,urot,theta)
            #####################################
            # check for symmetric molecule
            if distance(substr.getAtom(substreact).coords(),substr.centersym()) < 8.0e-2:
                atsc = substr.getBondedAtoms(substreact)
                r0a = substr.getAtom(substreact).coords()
                r1a = substr.getAtom(atsc[0]).coords()
                r2a = substr.getAtom(atsc[1]).coords()
                theta,u = rotation_params(r0a,r1a,r2a)
                theta = vecangle(u,vecdiff(r0a,mcoords))
                urot = cross(u,vecdiff(r0a,mcoords))
                ####################################
                # rotate around axis and get both images
                substrb = mol3D()
                substrb.copymol3D(substr)
                substr = rotate_around_axis(substr,r0a,urot,theta)
                substrb = rotate_around_axis(substrb,r0a,urot,-theta)
                d2 = substr.mindist(ts3D)
                d1 = substrb.mindist(ts3D)
                substr = substr if (d1 < d2)  else substrb # pick best one
            # rotate around axis of symmetry and get best orientation
            r1 = substr.getAtom(substreact).coords()
            u = vecdiff(r1,mcoords)
            dtheta = 2
            optmax = -9999
            totiters = 0
            substrb = mol3D()
            substrb.copymol3D(substr)
            # check for minimum distance between atoms and center of mass distance
            while totiters < 180:
                substr = rotate_around_axis(substr,r1,u,dtheta)
                d0 = substr.mindist(core) # try to maximize minimum atoms distance
                d0cm = substr.distance(core) # try to maximize center of mass distance
                iteropt = d0cm+10*log(d0) # optimization function
                if (iteropt > optmax): # if better conformation, keep
                    substrb = mol3D()
                    substrb.copymol3D(substr)
                    optmax = iteropt
                totiters += 1
            substr = substrb
    # combine molecules
    ts3D = mol3D()
    ts3D.copymol3D(core)
    ts3D = ts3D.combine(substr)
    ts3D.sanitycheck(0)
    ############ END FUNCTIONALIZING ###########
    fname = name_TS(rootdir,args.core,substr,args,bind=args.bind,bsmi=args.nambsmi)
    ts3D.writexyz(fname)
    strfiles.append(fname)
    getinputargs(args,fname)
    pfold = rootdir.split('/',1)[-1]
    # check for molecule sanity
    sanity,d0 = ts3D.sanitycheck(True)
    if args.debug:
        print('setting sanity diag, min dist at ' +str(d0) + ' (higher is better)')
    this_diag.set_sanity(sanity,d0)
    this_diag.set_mol(ts3D)
    this_diag.write_report(fname+'.report')
    del ts3D
    if sanity:
        print 'WARNING: Generated complex is not good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
    print '\nIn folder '+pfold+' generated 1 structure(s)!'
    return strfiles, emsg, this_diag
