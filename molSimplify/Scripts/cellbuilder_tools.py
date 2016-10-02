# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT
import os, sys, copy
import glob, re, math, random, string, numpy, pybel
from math import pi
from scipy.spatial import Delaunay, ConvexHull
#import networkx as nx
from molSimplify.Scripts.geometry import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.mol3D import*
from molSimplify.Classes.globalvars import globalvars
from operator import add
from molSimplify.Scripts.periodic_QE import *
###############################
def cell_ffopt(ff,mol,frozenats):
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
    ### convert mol3D to OBmol via xyz file, because AFTER/END option have coordinates
    backup_mol = mol3D()
    backup_mol.copymol3D(mol)
 #   print('bck ' + str(backup_mol.getAtom(0).coords()))
 #   print('mol_ibf ' + str(mol.getAtom(0).coords()))

    mol.writexyz('tmp.xyz')
    mol.OBmol = mol.getOBmol('tmp.xyz','xyzf')
    os.remove('tmp.xyz')
    ### initialize constraints
    constr = pybel.ob.OBFFConstraints()
    ### openbabel indexing starts at 1 ### !!!
    # convert metals to carbons for FF
    indmtls = []
    mtlsnums = []
    for iiat,atom in enumerate(mol.OBmol.atoms):
        if atom.atomicnum in metals:
            indmtls.append(iiat)
            mtlsnums.append(atom.atomicnum)
            atom.OBAtom.SetAtomicNum(19)
    for cat in frozenats:
        constr.AddAtomConstraint(cat+1) # indexing babel
    ### set up forcefield
    forcefield =pybel.ob.OBForceField.FindForceField(ff)
    obmol = mol.OBmol.OBMol
    forcefield.Setup(obmol,constr)
    ## force field optimize structure
    forcefield.ConjugateGradients(2500)
    forcefield.GetCoordinates(obmol)
    mol.OBmol = pybel.Molecule(obmol)

    # reset atomic number to metal
    for i,iiat in enumerate(indmtls):
        mol.OBmol.atoms[iiat].OBAtom.SetAtomicNum(mtlsnums[i])
    mol.convert2mol3D()

    en = forcefield.Energy()
 #   print(str(mol.OBmol.atoms[1].OBAtom.GetVector().GetZ()))
#    print(str(forcefield.Validate()))
   # print('mol_af ' + str(mol.getAtom(0).coords()))

  #  print('ff delta = ' + str(backup_mol.rmsd(mol)))
    del forcefield, constr, obmol
    return mol,en
################################
def import_from_cif(fst):
    #INPUT:
    # fst:  filename of cif file
    #OUTPUT:
    # unit_cell:  mol3D class of a single unit cell
    # cell_vector: list of lists of floats, each 
    #           corresponds to one of the defining cell
    #           vectors 
    cell_vector = list()
    unit_cell = mol3D()
    A = 0
    B = 0
    C = 0
    alpha =0
    beta = 0
    emsg =list()
    exit_status = 0
    gamma = 0
    obConversion = pybel.ob.OBConversion()
    obConversion.SetInAndOutFormats("cif", "xyz")
    mol = pybel.ob.OBMol()
    try:
        obConversion.ReadFile(mol, fst)
        fillUC = pybel.ob.OBOp.FindType("fillUC")
        fillUC = pybel.ob.OBOp.FindType("fillUC")
        fillUC.Do(mol, "strict")
        unit_cell.OBmol =pybel.Molecule(mol)
        unit_cell.convert2mol3D()
    except:
        emsg.append("Error in reading of cif file by pybel")
        exit_status = 1
    with open(fst) as f:
        lines = f.readlines()
        for line in lines:
            linesplit = line.split()
            if len(linesplit) != 0:
                if linesplit[0] == "_cell_length_a":
                    A = float((re.sub(r'\([^)]*\)','', ''.join(c for c in linesplit[1]))))
                if linesplit[0] == "_cell_length_b":
                   B = float((re.sub(r'\([^)]*\)','', ''.join(c for c in linesplit[1]))))
                if linesplit[0] == "_cell_length_c":
                    C = float((re.sub(r'\([^)]*\)','', ''.join(c for c in linesplit[1]))))

                if linesplit[0] == "_cell_angle_alpha":
                    alpha =float( ''.join(c for c in linesplit[1] if c not in '()').rstrip('.'))
                if linesplit[0] == "_cell_angle_beta":
                    beta =float( ''.join(c for c in linesplit[1] if c not in '()').rstrip('.'))
                if linesplit[0] == "_cell_angle_gamma":
                    gamma = float(''.join(c for c in linesplit[1] if c not in '()').rstrip('.'))
    # create cell vectors
    try:
        cell_vector.append([A,0,0])
        cell_vector.append([B*numpy.cos((gamma*pi)/180),B*numpy.sin((gamma*pi)/180),0])
        cx = C*numpy.cos((beta*pi)/180)
        cy = C*(numpy.cos((alpha*pi)/180)-numpy.cos((beta*pi)/180)*numpy.cos((gamma*pi/180)))/numpy.sin((gamma*pi)/180)
        cz = sqrt(C*C - cx*cx - cy*cy)
        cell_vector.append([cx,cy,cz])
    except:
        emsg = emsg.append('Error in creating unit cell from cif informtation')
        exit_status = 2
    for i,rows in enumerate(cell_vector):
        print(rows)
        for j,elements in enumerate(rows):
            if elements <= 1e-8:
                cell_vector[i][j] = 0
    if exit_status != 0:
        return emsg
    else:
        return unit_cell,cell_vector
##################################
def center_of_sym(list_of_points):
    n = len(list_of_points)
#    print('lop = ' + str(list_of_points))
    csym = [0,0,0];
    csym = [float(sum(x)/n) for x in zip(*list_of_points)]
    return csym
def zero_z_csm(super_cell):
    # puts center of sym
    # at z = 0
    csm = super_cell.centersym()
    vec = [0,0,0]
    vec[2] = -1*csm[2]
    super_cell.translate(vec)
##################################
def xgcd(b, n):
    # calculate x,y such that b*x + n*y = gcd(b,n)
    # by extended Euclidean algorithm 
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return  b, x0, y0
#
##################################

def distance_zw(r1,r2):
    dx = r1[0] - r2[0]
    dy = r1[1] - r2[1]
    dz =150*( r1[2] - r2[2])
    d = sqrt(dx**2+dy**2+dz**2)
    return d
##################################
def mdistance(r1,r2):
    dx = r1[0] - r2[0]
    dy = r1[1] - r2[1]
    dz = r1[2] - r2[2]
    d = sqrt(numpy.power(dx,2) + numpy.power(dy,2) + numpy.power(dz,2))
    return d

###################################
def normalize_vector(v):
    length = distance(v,[0,0,0])
    if length:
        nv = [float(i)/length for i in v]
    else:
        nv = [0,0,0]
    return nv
##############################
def find_all_surface_atoms(super_cell,tol=1e-2,type_of_atom = False):
    # Get all atoms on the tope surface - NB, this will
    # not handle complex (2 or more) atom-type surfaces
    # if the atoms are 'layered', e.g. TiO2 - Ti under O2, 
    # no Ti will be found. This can be overcome by using \
    # a looser tolerance, such that the interlayer differences
    # are smaller than tol, but not so large as to conflate 
    # different layers!
    # INPUT: 
    #   - super_cell: mol3D class that contains the super cell
    #   - tol: float, max distance from extent plane to look
    #   - type_of_atom: optional, string, gets atoms of the given type on the face plane
    #                   if left out, will not care about types of atoms
    # OUPUT
    #   - avail_sites_list: list of int, indices of atoms on the surface
    #
    extents = find_extents(super_cell)
    target_height = extents[2]
    avail_sites_list = list()
    if type_of_atom:
        possible_atom_inds = super_cell.findAtomsbySymbol(type_of_atom)
    else:
        possible_atom_inds =range(0,super_cell.natoms)
    for indices in possible_atom_inds:
        z_dist = abs(super_cell.getAtom(indices).coords()[2] - target_height)
        if (z_dist <= tol):
            avail_sites_list.append(indices)

    return avail_sites_list

###################################
def distance_2d_torus(R1,R2,dim):
    ### distance between points in Euclidean torus 
    dx =abs( R1[0] - R2[0])
    dy = abs(R1[1] - R2[1] )
    dz =abs(( R1[2] - R2[2]))
    d1 = sqrt(  numpy.power(dim[0] - dx,2)
              + numpy.power(dim[1] - dy,2)
              + numpy.power(dz,2))
    d2 = sqrt(  numpy.power(dim[0] - dx,2)
                     + numpy.power(dy,2)
                     + numpy.power(dz,2))
    d3 = sqrt(  numpy.power(dx,2)
              + numpy.power(dim[1] - dy,2)
              + numpy.power(dz,2))
    d = min(d1,d2,d3)
    return d
################################################################
def periodic_2d_distance(R1,R2,cell_vector):
    ### distance between points in Euclidean torus 
    ## STILL UNDER CONSTRUCTION, WIP WIP WIP ***
    dx =abs( R1[0] - R2[0])
    dy = abs(R1[1] - R2[1] )
    dz =abs(( R1[2] - R2[2]))
    for v1shifts in [-1,0,1]:
        for v1shifts in [-1,0,1]:
            for yshifts in [-1,0,0]:
                pass
    d1 = sqrt(  numpy.power(dim[0] - dx,2)
              + numpy.power(dim[1] - dy,2)
              + numpy.power(dz,2))
    d2 = sqrt(  numpy.power(dim[0] - dx,2)
                     + numpy.power(dy,2)
                     + numpy.power(dz,2))
    d3 = sqrt(  numpy.power(dx,2)
              + numpy.power(dim[1] - dy,2)
              + numpy.power(dz,2))
    d = min(d1,d2,d3)
    return d
################################################################
def periodic_mindist(mol,surf,dim):
    ### calculates minimum distance between atoms in 2 molecules ###
    # INPUT
    #   - mol: mol3D class,  molecule
    #   - surf: mol3D class, the surface
    #   - dim: list of float, replication 
    # OUTPUT
    #   - mind: minimum distance between atoms of the 2 mol objects
    mind = 1000
    for atom1 in mol.getAtoms():
        for atom0 in surf.getAtoms():
            if (distance_2d_torus(atom1.coords(),atom0.coords(),dim) < mind):
                mind = distance(atom1.coords(),atom0.coords())
    return mind
################################################################
def periodic_selfdist(mol,dim):
    ### calculates minimum distance between atoms in 2 molecules ##
    # INPUT
    #   - mol: mol3D class,  molecule
    #   - dim: list of floats, replication 
    # OUTPUT
    #   - mind: minimum distance between atoms of the 2 mol and periodic
    #             images
    mind = 1000
    for ii,atom1 in enumerate(mol.getAtoms()):
        for jj,atom0 in enumerate(mol.getAtoms()):
            if (distance_2d_torus(atom1.coords(),atom0.coords(),dim) < mind) and (ii !=jj):
                mind = distance(atom1.coords(),atom0.coords())
    return mind

##################################
def closest_torus_point(mol,dim):
    min_dist = 1000
    for atom1 in mol.getAtoms():
        R1 = atom1.coords()
        for atom2 in mol.getAtoms():
            R2 = atom2.coords()
            d = distance_2d_torus(R1,R2,dim)
            if (d<min_dist):
                min_dist = d
    return min_dist
#################################
def check_top_layer_correct(super_cell,atom_type):
    # remove the layer on 
    # top of the cell if
    # wrong material is exposed
    trimmed_cell = mol3D()
    trimmed_cell.copymol3D(super_cell)
    globs = globalvars()
    elements = globs.elementsbynum()
    print('chekcing surface  for  ' + atom_type + '\n')
    if not atom_type in elements:
        print("unkown surface type, unable to trim ")
        return trimmed_cell
    else:
        stop_flag = False
        counter = 0 # 3 tries max
        while not stop_flag:
            atom_type_surf = find_all_surface_atoms(trimmed_cell,tol=1e-3,type_of_atom = atom_type)
            top_surf = find_all_surface_atoms(trimmed_cell,tol=1e-3,type_of_atom = False)
#            print("top surf",top_surf)
#            print("atom top surf",atom_type_surf)
            if set(atom_type_surf) == set(top_surf):
                print('match')
                stop_flag = True
            else:
                counter += 1
                trimmed_cell = shave_surface_layer(trimmed_cell,1e-3)
            if counter == 3:
                print('unable to find target atom in 3 cuts')
                stop_flag = True
    return trimmed_cell
###############################

def shave_surface_layer(super_cell,TOL=1e-1):
    shaved_cell = mol3D()
    shaved_cell.copymol3D(super_cell)
    extents = find_extents(super_cell)
    zmax = extents[2]
    del_list = list()
    for i,atoms in enumerate(super_cell.getAtoms()):
        coords = atoms.coords()
        if abs(coords[2] - zmax) < TOL:
            del_list.append(i)
    shaved_cell.deleteatoms(del_list)
    return shaved_cell
###############################
def shave_under_layer(super_cell):
    shaved_cell = mol3D()
    shaved_cell.copymol3D(super_cell)
    TOL = 1e-1
    zmin = 1000;
    for i,atoms in enumerate(super_cell.getAtoms()):
        coords = atoms.coords()
        if (coords[2] < zmin):
                zmin = coords[2]
    del_list = list()
    for i,atoms in enumerate(super_cell.getAtoms()):
        coords = atoms.coords()
        if abs(coords[2] - zmin) < TOL:
            del_list.append(i)
    shaved_cell.deleteatoms(del_list)
    return shaved_cell
###############################
def shave_dim_type(super_cell,dim,mode):
    ## dim  = 0,1,2
    ##        x,y,z 
    #   mode = 1 for max, -1 for min
    shaved_cell = mol3D()
    shaved_cell.copymol3D(super_cell)
    TOL = 1e-1

    dim_ref = 1000

    if (mode == -1):
        for i,atoms in enumerate(super_cell.getAtoms()):
            coords = atoms.coords()
            if (coords[dim] < dim_ref):
                    dim_ref = coords[dim]
        del_list = list()
        for i,atoms in enumerate(super_cell.getAtoms()):
            coords = atoms.coords()
            if abs(coords[dim] - dim_ref) < TOL:
                del_list.append(i)
    if (mode ==1):
        extents = find_extents(super_cell)
        dim_max = extents[dim]
        del_list = list()
        for i,atoms in enumerate(super_cell.getAtoms()):
            coords = atoms.coords()
            if abs(coords[dim] - dim_max) < TOL:
                del_list.append(i)
    ## return
    shaved_cell.deleteatoms(del_list)
    return shaved_cell
###############################

def zero_z(super_cell):
        zeroed_cell = mol3D()
        zeroed_cell.copymol3D(super_cell)
        TOL = 1e-1
        zmin = 1000;
        for i,atoms in enumerate(super_cell.getAtoms()):
                coords = atoms.coords()
                if (coords[2] < zmin):
                        zmin = coords[2]
        zeroed_cell.translate([0,0,-1*zmin])
        return zeroed_cell
def zero_x(super_cell):
        zeroed_cell = mol3D()
        zeroed_cell.copymol3D(super_cell)
        TOL = 1e-1
        xmin = 1000;
        for i,atoms in enumerate(super_cell.getAtoms()):
                coords = atoms.coords()
                if (coords[0] < xmin):
                        xmin = coords[0]
        zeroed_cell.translate([-1*xmin,0,0])
        return zeroed_cell
def zero_y(super_cell):
        zeroed_cell = mol3D()
        zeroed_cell.copymol3D(super_cell)
        TOL = 1e-1
        ymin = 1000;
        for i,atoms in enumerate(super_cell.getAtoms()):
                coords = atoms.coords()
                if (coords[1] < ymin):
                        ymin = coords[1]
        zeroed_cell.translate([0,-1*ymin,0])
        return zeroed_cell

###############################
def point_in_box(point,box):
    outcome = False
    fx=(box[0][0] <= point[0])*(point[0] < box[0][1])
    fy=(box[1][0] <= point[1])*(point[1] < box[1][1])
    fz=(box[2][0] <= point[2])*(point[2] < box[2][1])
    if fz  and fy and fx:
        outcome = True
    return outcome
##############################
def points_below_plane(point,w):
    outcome = False
    zplane = w[2] + w[1]*point[1] + w[0]*point[0]
    if (point[2] <= zplane):
        outcome =  True
    return outcome
#################################
def freeze_bottom_n_layers(super_cell,n):
        frozen_cell = mol3D()
        frozen_cell.copymol3D(super_cell)
        counter = 0 
        while  counter < n:
                frozen_cell_new = freeze_under_layer(frozen_cell)
                frozen_cell = mol3D()
                frozen_cell.copymol3D(frozen_cell_new)
                counter +=1
        return frozen_cell
###############################
def freeze_under_layer(super_cell):
    frozen_cell = mol3D()
    frozen_cell.copymol3D(super_cell)
    TOL = 0.5
    zmin = 1000;
    for i,atoms in enumerate(super_cell.getAtoms()):
        coords = atoms.coords()
        if not atoms.frozen:
                if (coords[2] < zmin):
                        zmin = coords[2]
    freeze_list = list()
 #   print('lowest  ' + str(zmin))
    for i,atoms in enumerate(super_cell.getAtoms()):
        coords = atoms.coords()
        if not atoms.frozen:
                if abs(coords[2] - zmin) < TOL:
                    freeze_list.append(i)
#                    print("freezing at " + str(coords))
    frozen_cell.freezeatoms(freeze_list)
    return frozen_cell
##############################
def find_extents(super_cell):
    # INPUT
    #   - super_cell: mol3D class that contains the super cell
    # OUPUT
    #   - extents: list of max coords of atoms on the surface
    xmax = 0
    zmax = 0
    ymax = 0
    for atoms in super_cell.getAtoms():
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
def find_extents_cv(super_cell_vector):
    # INPUT
    #   - super_cell_vector: matrix of the three vectors that define the super cell
    # OUPUT
    #   - extents: list of max coords of the super cell
    xmax = 0
    zmax = 0
    ymax = 0
    for columns in super_cell_vector:
        xmax = max(xmax,abs(columns[0]))
        ymax = max(ymax,abs(columns[1]))
        zmax = max(zmax,abs(columns[2]))
    xmax = numpy.linalg.norm(super_cell_vector[0])
    ymax = numpy.linalg.norm(super_cell_vector[1])
    zmax = numpy.linalg.norm(super_cell_vector[2])

    extents = [xmax,ymax,zmax]
    return extents

