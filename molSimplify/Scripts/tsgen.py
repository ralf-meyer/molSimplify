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
#from pybel import *
#import openbabel
from pkg_resources import resource_filename, Requirement
import random, itertools
from numpy import log, arccos, cross, dot, pi

MOthrc = 0.9 # threshold for detection of M=O bond as a multiple of sum cov rad
#TS_HX_dist = 1.8 # TS geometric params in angstroms (to be trained)
#TS_OH_dist = 1.0
#TS_MO_dist = 1.76
#angopt = 140
TS_HX_dist = 1.8
TS_OH_dist = 0.98
TS_MO_dist = 1.76
angopt = 160

def getconnection(core,catom,Midx,BL,angopt):
    angcoeff = 0.1 # coefficient of MOH term in connecting point optimization
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

def getconnection2(core,Midx,BL):
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
			# objective function is based on maximizing the smallest angle with other connecting atoms
			d = []
			for iligs in nligs:
				d.append(vecangle(vecdiff(Mcoords,P),vecdiff(Mcoords,core.getAtomCoords(iligs))))
			d0 = min(d)
			if d0 > mdist:
				mdist = d0
				cpoint = P
	return cpoint

def checkMO(core,MOthr,Midx):
	MOcheck = 0 # count number of M=O bonds
	for i in core.getBondedAtoms(Midx):
		if core.getAtom(i).symbol() == 'O' and core.getAtom(i).distance(core.getAtom(Midx)) < MOthr:
			Oidx = i
			MOcheck += 1
	if MOcheck == 0:
		print('Error: No M=O bond could be found.\n')
	elif MOcheck > 1:
		print('Warning: Multiple M=O bonds detected. Using the last one.\n')
	elif MOcheck == 1:
		print('Generating TS guess...\n')
	return Oidx



def tsgen1(core,subst,nH):
	# load core
	ts3D = mol3D()
	ts3D.readfromxyz(core)
	Midx = ts3D.findMetal()
	MOthr = MOthrc*(ts3D.getAtom(Midx).rad + atom3D('O').rad) # M=O bond detection threshold
	Oidx = checkMO(ts3D,MOthr,Midx)
	# stretch M-O bond according to trained params
	MOnew = TS_MO_dist
	ts3D.BCM(Oidx,Midx,MOnew)
	# set O-H bond distance from trained params
	OHnew = TS_OH_dist
	Ocoords = ts3D.getAtom(Oidx).coords()
	# get connecting point of abstracted atom
	cpoint = getconnection(ts3D,Oidx,Midx,OHnew,angopt)
	ts3D.addAtom(atom3D('X',cpoint))
	#ts3D.writexyz('test.xyz')
	ts3D.deleteatom(ts3D.natoms-1)
	# load substrate
	lig3D = mol3D()
	lig3D.readfromxyz(subst)
	# stretch substrate H-X bond according to trained params
	adjligatidx = lig3D.getBondedAtoms(nH)
	if len(adjligatidx) > 1:
		print('Warning: Abstracted H is not terminal!\n')
	elif len(adjligatidx) == 0:
		print('Warning: Abstracted H is not bonded to anything!\n')
	else:
		(adjligatidx,) = adjligatidx
	HXnew = TS_HX_dist
	at = lig3D.getBondedAtoms(nH)
	lig3D.BCM(nH,adjligatidx,HXnew)
	# align substrate according to connection atom and shadow atom
	lig3D.alignmol(lig3D.getAtom(nH),atom3D('H',cpoint))
	# perform rotations
	mcoords = ts3D.getAtomCoords(Oidx)
	r0 = ts3D.getAtomCoords(Oidx)
	r1 = lig3D.getAtomCoords(nH)
	r2 = lig3D.centermass() # center of mass
	if not r2:
		emsg = 'Center of mass calculation for substrate failed. Check input.'
		print emsg
	rrot = r1
	theta,u = rotation_params(r0,r1,r2)
	# for most ligands align center of mass of local environment
	if (lig3D.natoms > 1):
		lig3Db = mol3D()
		lig3Db.copymol3D(lig3D)
		####################################
		# center of mass of local environment (to avoid bad placement of bulky ligands)
		auxmol = mol3D()
		for i in at:
			auxmol.addAtom(lig3D.getAtom(i))
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
		for i in at: # BUGGY, SHOULD BE CHECKED IN THE PRE-DISTORTED SUBSTRATE INSTEAD
			auxm.addAtom(lig3D.getAtom(i))
		if auxm.natoms > 1:
			r0 = lig3D.getAtom(nH).coords()
			r1 = auxm.getAtom(0).coords()
			r2 = auxm.getAtom(1).coords()
			if checkcolinear(r1,r0,r2):
				# we will rotate so that angle is right
				theta,urot = rotation_params(r1,mcoords,r2)
				theta = vecangle(vecdiff(r0,mcoords),urot)
				lig3D = rotate_around_axis(lig3D,r0,urot,theta)
		#####################################
		# check for symmetric molecule
		if distance(lig3D.getAtom(nH).coords(),lig3D.centersym()) < 8.0e-2:
			atsc = lig3D.getBondedAtoms(nH)
			r0a = lig3D.getAtom(nH).coords()
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
			d2 = lig3D.mindist(ts3D)
			d1 = lig3Db.mindist(ts3D)
			lig3D = lig3D if (d1 < d2)  else lig3Db # pick best one
		# rotate around axis of symmetry and get best orientation
		r1 = lig3D.getAtom(nH).coords()
		u = vecdiff(r1,mcoords)
		dtheta = 2
		optmax = -9999
		totiters = 0
		lig3Db = mol3D()
		lig3Db.copymol3D(lig3D)
		# check for minimum distance between atoms and center of mass distance
		while totiters < 180:
			lig3D = rotate_around_axis(lig3D,r1,u,dtheta)
			d0 = lig3D.mindist(ts3D) # try to maximize minimum atoms distance
			d0cm = lig3D.distance(ts3D) # try to maximize center of mass distance
			iteropt = d0cm+10*log(d0) # optimization function
			if (iteropt > optmax): # if better conformation, keep
				lig3Db = mol3D()
				lig3Db.copymol3D(lig3D)
				optmax = iteropt
			totiters += 1
		lig3D = lig3Db
	# get correct distance for center of mass
	#u = vecdiff(cpoint,Ocoords)
	#lig3D = aligntoaxis2(lig3D, cpoint, Ocoords, u, bondl)
	# combine molecules
	ts3D.deleteatom(ts3D.natoms-1)
	ts3D = ts3D.combine(lig3D)
	ts3D.writexyz('TS.xyz')
	ts3D.sanitycheck(0)

def tsgen2(core,subst,nH):
	MXYopt = 95
	# load core and substrate
	ts3D = mol3D()
	ts3D.readfromxyz(core)
	Midx = ts3D.findMetal()
	# check for empty coordination site
	if len(ts3D.getBondedAtoms(Midx)) > 5:
		print('Error, no empty coordination site found, exiting...')
		sys.exit()
	lig3D = mol3D()
	lig3D.readfromxyz(subst)
	# stretch substrate bond according to trained params
	adjligatidx = lig3D.getBondedAtoms(nH)
	if len(adjligatidx) > 1:
		print('Warning: Bonding atom is not terminal!\n')
	elif len(adjligatidx) == 0:
		print('Warning: Bonding atom is not bonded to anything!\n')
	else:
		(adjligatidx,) = adjligatidx
	XYnew = TS_XY_dist*(lig3D.getAtom(nH).rad + lig3D.getAtom(adjligatidx).rad)
	lig3D.BCM(nH,adjligatidx,XYnew)
	# place substrate at empty coordination site
	BL = TS_MX_dist*(ts3D.getAtom(Midx).rad + lig3D.getAtom(nH).rad)
	cpoint = getconnection2(ts3D,Midx,BL)
	lig3D.alignmol(lig3D.getAtom(nH),atom3D('H',cpoint))
	# align substrate to desired MXY angle
	tmp3D = mol3D()
	tmp3D.copymol3D(ts3D)
	tmp3D.addAtom(atom3D('Cl',cpoint))
	ligalignpt = getconnection(tmp3D,tmp3D.natoms-1,Midx,BL,MXYopt)
	theta,u = rotation_params(lig3D.getAtomCoords(adjligatidx),cpoint,ligalignpt)
	lig3D = rotate_around_axis(lig3D,cpoint,u,180-theta)
	# for ligands with more than 2 atoms, rotate around X-Y axis to minimize steric repulsion
	if lig3D.natoms > 2:
		dtheta = 2
		optmax = -9999
		totiters = 0
		lig3Db = mol3D()
		lig3Db.copymol3D(lig3D)
		# check for minimum distance between atoms and center of mass distance
		while totiters < 180:
			lig3D = rotate_around_axis(lig3D,cpoint,vecdiff(cpoint,ligalignpt),dtheta)
			d0 = lig3D.mindist(ts3D) # try to maximize minimum atoms distance
			d0cm = lig3D.distance(ts3D) # try to maximize center of mass distance
			iteropt = d0cm+10*log(d0) # optimization function
			if (iteropt > optmax): # if better conformation, keep
				lig3Db = mol3D()
				lig3Db.copymol3D(lig3D)
				optmax = iteropt
			totiters += 1
		lig3D = lig3Db
	# combine molecules
	ts3D = ts3D.combine(lig3D)
	ts3D.writexyz('TS.xyz')
	ts3D.sanitycheck(0)		

def tsgen(args,rootdir,substr,globs):
    # INPUT
    #   - args: placeholder for input arguments
    #   - rootdir: directory of current run
    #   - substr: substrate mol3D
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
    # load molecule core
    core,emsg = core_load(args.core,mcores)
    if emsg:
        return False,emsg
    core.convert2mol3D() # convert to mol3D
    # copy initial core for backup
    initcore3D = mol3D()
    initcore3D.copymol3D(core)
    sanity = False
    # placeholder
    core3D = initcore3D
    ############ END FUNCTIONALIZING ###########
    fname = name_TS(rootdir,args.core,substr,args,bind=args.bind,bsmi=args.nambsmi)
    core3D.writexyz(fname)
    strfiles.append(fname)
    getinputargs(args,fname)
    pfold = rootdir.split('/',1)[-1]
    # check for molecule sanity
    sanity,d0 = core3D.sanitycheck(True)
    if args.debug:
        print('setting sanity diag, min dist at ' +str(d0) + ' (higher is better)')
    this_diag.set_sanity(sanity,d0)
    this_diag.set_mol(core3D)
    this_diag.write_report(fname+'.report')
    del core3D
    if sanity:
        print 'WARNING: Generated complex is not good! Minimum distance between atoms:'+"{0:.2f}".format(d0)+'A\n'
    print '\nIn folder '+pfold+' generated 1 structure(s)!'
    return strfiles, emsg, this_diag
