
# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

import sys, os, random, shutil, inspect, argparse, openbabel
from molSimplify.Scripts.rungen import *
from molSimplify.Scripts.io import *
from molSimplify.Scripts.inparse import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.mol3D import*
from molSimplify.Classes.globalvars import globalvars
from math import sqrt
from math import floor
def chain_ffopt(ff,mol,frozenats):
		### FORCE FIELD OPTIMIZATION ##
		# INPUT
		#   - ff: force field to use, available MMFF94, UFF< Ghemical, GAFF
		#   - mol: mol3D to be ff optimized
		#   - connected: indices of connection atoms to metal
		#   - constopt: flag for constrained optimization
		# OUTPUT
		#   - mol: force field optimized mol3D
		metals = range(21,31)+range(39,49)+range(72,81)
		### check requested force field
		ffav = 'mmff94, uff, ghemical, gaff, mmff94s' # force fields
		if ff.lower() not in ffav:
			print 'Requested force field not available. Defaulting to MMFF94'
			ff = 'mmff94'
		### convert mol3D to OBMol via xyz file, because AFTER/END option have coordinates
		backup_mol = mol3D()
		backup_mol.copymol3D(mol)
		#   print('bck ' + str(backup_mol.getAtom(0).coords()))
		#   print('mol_ibf ' + str(mol.getAtom(0).coords()))
		mol.convert2OBMol()
		### initialize constraints
		constr = openbabel.OBFFConstraints()
		### openbabel indexing starts at 1 ### !!!
		# convert metals to carbons for FF
		indmtls = []
		mtlsnums = []
		for iiat,atom in enumerate(openbabel.OBMolAtomIter(OBMol)):
			if atom.atomicnum in metals:
				indmtls.append(iiat)
				mtlsnums.append(atom.GetAtomicNum())
				atom.OBAtom.SetAtomicNum(19)
		for cat in frozenats:
			constr.AddAtomConstraint(cat+1) # indexing babel
		### set up forcefield
		forcefield = openbabel.OBForceField.FindForceField(ff)
		OBMol = mol.OBMol
		forcefield.Setup(obmol,constr)
		## force field optimize structure
		forcefield.ConjugateGradients(2500)
		forcefield.GetCoordinates(OBMol)
		mol.OBMol = OBMol

# reset atomic number to metal
		for i,iiat in enumerate(indmtls):
			mol.OBMol.GetAtomById(iiat).SetAtomicNum(mtlsnums[i])
		mol.convert2mol3D()

		en = forcefield.Energy()
		#   print(str(mol.OBmol.atoms[1].OBAtom.GetVector().GetZ()))
		#    print(str(forcefield.Validate()))
		# print('mol_af ' + str(mol.getAtom(0).coords()))

		#  print('ff delta = ' + str(backup_mol.rmsd(mol)))
		del forcefield, constr, OBMol
		return mol,en

def mdistance(r1,r2):
	dx = r1[0] - r2[0]
	dy = r1[1] - r2[1]
	dz = r1[2] - r2[2]
	d = sqrt(numpy.power(dx,2) + numpy.power(dy,2) + numpy.power(dz,2))
	return d
def find_extents(mol):
	# INPUT
	#   - mol: mol3D class that contains the super cell
	# OUPUT
	#   - extents: list of max coords of atoms on the surface
	xmax = 0
	zmax = 0
	ymax = 0
	for atoms in mol.getAtoms():
		coords = atoms.coords()
		x_ext = coords[0]# + atoms.rad
		y_ext = coords[1]# + atoms.rad
		z_ext = coords[2]# + atoms.rad
		xmax = max(xmax,x_ext)
		ymax = max(ymax,y_ext)
		zmax = max(zmax,z_ext)
	extents = [xmax,ymax,zmax]
	return extents
#####################################
def interatomic_dist(mol,ind1,ind2):
	   coord1 = mol.getAtom(ind1).coords()
	   coord2 = mol.getAtom(ind2).coords()
	   distance =  mdistance(coord1,coord2)
	   print('iad between',mol.getAtom(ind1).symbol(),mol.getAtom(ind2).symbol())

	   vector  =  [coord1[i] - coord2[i] for i in [0,1,2]]
	   return distance,vector
def find_term_heavy(mol,reference_point):
		min_dist  = 1000
		min_ind = 0
	#    print('reference_point',reference_point)
		for inds,atoms in enumerate(mol.getAtoms()):
	 #           print('inds,atoms')
	  #          print('this ind is ',inds)
	  #          print('this atom ',atoms.symbol())
				if not atoms.symbol() == "H":
					   this_coord = atoms.coords()
					   this_dist = distance(this_coord,reference_point)
	   #                print(atoms.symbol(),inds,this_coord,this_dist)
					   if this_dist < min_dist:
							   min_dist = this_dist
							   min_ind = inds
							   print('accepting')
		return min_ind
def trim_H(mol,reference_point):
		trimmed_mol = mol3D()
		trimmed_mol.copymol3D(mol)
		min_ind = find_term_heavy(mol,reference_point)
		hydrogen_list = trimmed_mol.getHsbyIndex(min_ind)
		trimmed_mol.deleteatoms([hydrogen_list[0]])
		return trimmed_mol
def zero_dim(mol,dim):
		if dim == 0:
				return zero_x(mol)
		elif dim == 1:
				return zero_y(mol)
		elif dim == 2:
				return zero_z(mol)
		else:
				return 'Error'
def zero_z(mol):
		zeroed_mol = mol3D()
		zeroed_mol.copymol3D(mol)
		TOL = 1e-1
		zmin = 1000;
		for i,atoms in enumerate(mol.getAtoms()):
				coords = atoms.coords()
				if (coords[2] < zmin):
						zmin = coords[2]
		zeroed_mol.translate([0,0,-1*zmin])
		return zeroed_mol
def zero_x(mol):
		zeroed_mol = mol3D()
		zeroed_mol.copymol3D(mol)
		TOL = 1e-1
		xmin = 1000;
		for i,atoms in enumerate(mol.getAtoms()):
				coords = atoms.coords()
				if (coords[0] < xmin):
						xmin = coords[0]
		zeroed_mol.translate([-1*xmin,0,0])
		return zeroed_mol
def zero_y(mol):
		zeroed_mol = mol3D()
		zeroed_mol.copymol3D(mol)
		TOL = 1e-1
		ymin = 1000;
		for i,atoms in enumerate(mol.getAtoms()):
				coords = atoms.coords()
				if (coords[1] < ymin):
						ymin = coords[1]
		zeroed_mol.translate([0,-1*ymin,0])
		return zeroed_mol
def remove_closest_h(mol,other_mol):
	new_mol = mol3D()
	new_mol.copymol3D(mol);
	min_distance  = 1000
	current_ind = 0
	for Hatoms in mol.getHs():
		this_H = mol.getAtom(Hatoms)
		for atoms in other_mol.getAtoms():
			this_distance  = mdistance(this_H.coords(),atoms.coords())
			if this_distance < min_distance:
				min_distance = this_distance
				current_ind = Hatoms
	print(' the H ind to delete is  ' +str(current_ind) + '  at ' + str(min_distance)) 
	new_mol.deleteatoms([current_ind])
	return new_mol
	

def grow_linear_step(chain,new_unit,dim,interv):
		combined_mol = mol3D()
		combined_mol.copymol3D(chain)
		add_mol = mol3D()
		add_mol.copymol3D(new_unit)
		add_mol = zero_dim(new_unit,dim)
		chain_inds = range(0,chain.natoms)
		print('chain_inds',chain_inds)
		basic_lengths = find_extents(chain)
	   # print(basic_lengths)
		basic_dist  = basic_lengths[dim]
		tv =[0,0,0]
		tv[dim] = basic_dist
#        add_mol.translate(tv)
		print('translating',tv)

		add_mol.translate(interv)
		
                add_mol.writexyz('precut.xyz')
                add_mol = remove_closest_h(add_mol,combined_mol)
                add_mol.writexyz('postcut.xyz')
		combined_mol = combined_mol.combine(add_mol)
		#combined_mol,en = chain_ffopt('',combined_mol,[])       

		combined_mol.writexyz('pre.xyz')

		combined_mol.writexyz('post.xyz')

		return combined_mol

def chain_builder_supervisor(args,rundir):
		emsg = list()

		if not (args.chain) and not (isinstance(args.chain_units, (int))):
				emsg.append('Invalid input: need monomer AND number of units')

		print(args.chain)
		print(args.chain_units)

		monomer = mol3D()
		monomer.OBmol = monomer.getOBmol(args.chain,'smi')
		monomer.OBmol.make3D('mmff94',0)
		monomer.convert2mol3D()

		monomer.writexyz('mono.xyz')
		dimer = mol3D()
		dimer.OBmol = dimer.getOBmol(args.chain+args.chain,'smi')
		dimer.OBmol.make3D('mmff94',0)
		dimer.convert2mol3D()
		dimer.printxyz()
		dimer.writexyz('di.xyz')
                interd,interv = interatomic_dist(dimer,len(args.chain),0)
		print('interv is')
		print(interv)
                
		trimer = mol3D()
		trimer.OBmol = trimer.getOBmol(args.chain+args.chain + args.chain,'smi')
		trimer.OBmol.make3D('mmff94',0)
		trimer.convert2mol3D()
		trimer.printxyz()
		trimer.writexyz('tri.xyz')

		my_dim = mol3D()
		my_dim.copymol3D(monomer)
		my_dim.writexyz('prestart.xyz')
                
		my_dim = trim_H(my_dim,monomer.getAtom(len(args.chain)-1).coords())
		my_dim = zero_x(my_dim)
		basic_lengths = find_extents(my_dim)

		basic_x  = basic_lengths[0]
		basic_y  = basic_lengths[1]


		my_dim.writexyz('start.xyz')

		middle = mol3D()
		middle.copymol3D(monomer)
                print('connection atom is '+monomer.getAtom(len(args.chain)-1).symbol())
		middle = trim_H(middle,monomer.getAtom(len(args.chain)-1).coords())

		middle =zero_x(middle)
		middle = trim_H(middle,[0,0,0])

		end = mol3D()
		end.copymol3D(monomer)
		end =zero_x(end)
		end = trim_H(end,[0,0,0])

		middle.writexyz('middle.xyz')
		end.writexyz('end.xyz')



		repu = mol3D()
		repu.copymol3D(my_dim)







		interv0 = interv
		for i in range(0,int(args.chain_units)-1):
				my_dim = grow_linear_step(my_dim,repu,0,interv)
				interv = [interv[i] + interv0[i] for i in [0,1,2]]
		my_dim = grow_linear_step(my_dim,end,0,interv)

#        my_dim.printxyz()
		my_dim.writexyz('poly.xyz')
		my_dim,en = chain_ffopt('',my_dim,[])       
		my_dim.writexyz('polyf.xyz')


		if emsg:
				print(emsg)

		return emsg


