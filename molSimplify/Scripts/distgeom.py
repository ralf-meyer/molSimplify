## @file distgeom.py
#  Implements a basic distance geometry conformer search routine
#  
#  Written by Terry Gani for HJK Group
#  
#  Dpt of Chemical Engineering, MIT
# 
#  Adapted from: 
#  
#  [1] J. M. Blaney and J. S. Dixon, "Distance Geometry in Molecular Modeling", in Reviews in Computational Chemistry, VCH (1994)
#  
#  [2] G. Crippen and T. F. Havel, "Distance Geometry and Molecular Conformation", in Chemometrics Research Studies Series, Wiley (1988)

from molSimplify.Scripts.geometry import *
from molSimplify.Scripts.io import *
from molSimplify.Classes.globalvars import *
import os, sys
import openbabel
import numpy as np
from scipy import optimize
# for testing
import time

## Applies the cosine rule to get the length of AC given lengths of AB, BC and angle ABC
#  @param AB Length of AB
#  @param BC Length of BC
#  @param theta Angle in degrees
#  @return Length of AC
def CosRule(AB,BC,theta):
    theta = np.pi*theta/180
    AC = sqrt(AB**2+BC**2-2*AB*BC*cos(theta))
    return AC

## Generate distance bounds matrices
#  
#  The basic idea is outlined in ref [1].
#  
#  We first apply 1-2 (bond length) and 1-3 (bond angle) constraints, read from the FF-optimized initial conformer.
#  
#  Next, to bias the search towards coordinating conformers, approximate connection atom distance constraints based on topological distances are also included.
#  @param mol mol3D of molecule
#  @param natoms Number of atoms in molecule
#  @param catoms List of ligand connection atoms (default empty)
#  @param A Distance 2 connectivity matrix
#  @return Lower and upper bounds matrices
def GetBoundsMatrices(mol,natoms,catoms=[],A=[]):
    LB = np.zeros((natoms,natoms)) # lower bound
    UB = np.zeros((natoms,natoms)) # upper bound, both symmetric
    for i in range(natoms):
        for j in range(natoms):
            # 1-2 constraints: UB = LB = BL
            if mol.OBMol.GetBond(i+1,j+1) is not None:
                UB[i][j] = distance(mol.getAtomCoords(i),mol.getAtomCoords(j))
                UB[j][i] = distance(mol.getAtomCoords(i),mol.getAtomCoords(j))
                LB[i][j] = distance(mol.getAtomCoords(i),mol.getAtomCoords(j))
                LB[j][i] = distance(mol.getAtomCoords(i),mol.getAtomCoords(j))            
    for i in range(natoms):
        for j in range(natoms):
            for k in range(natoms):
            # 1-3 constraints: UB = LB = BL
                if mol.OBMol.GetBond(i+1,j+1) is not None and mol.OBMol.GetBond(j+1,k+1) is not None and j != k and i != k:
                    AB = vecdiff(mol.getAtomCoords(j),mol.getAtomCoords(i))
                    BC = vecdiff(mol.getAtomCoords(k),mol.getAtomCoords(j))
                    UB[i][k] = CosRule(norm(AB),norm(BC),180-vecangle(AB,BC))
                    UB[k][i] = CosRule(norm(AB),norm(BC),180-vecangle(AB,BC))      
                    LB[i][k] = CosRule(norm(AB),norm(BC),180-vecangle(AB,BC))
                    LB[k][i] = CosRule(norm(AB),norm(BC),180-vecangle(AB,BC))
    # set cis 1-2 constraints for connecting atoms
    for i,catom in enumerate(catoms[:-1]):
        # exclude catoms with distance 2 or less e.g. NCN
        if A[catoms[i],catoms[i+1]] == 0:
	        UB[catoms[i]][catoms[i+1]] = 2.8
	        UB[catoms[i+1]][catoms[i]] = 2.8
	        LB[catoms[i]][catoms[i+1]] = 2.8
	        LB[catoms[i+1]][catoms[i]] = 2.8
    # set right-triangle 1-3 constraints for connecting atoms
    if len(catoms) > 2:
        for i,catom in enumerate(catoms[:-2]):
	        UB[catoms[i]][catoms[i+2]] = 4.0
	        UB[catoms[i+2]][catoms[i]] = 4.0        
	        LB[catoms[i]][catoms[i+2]] = 4.0
	        LB[catoms[i+2]][catoms[i]] = 4.0
	# set planar 1-4 constraints for connecting atoms
    if len(catoms) > 3:
        for i,catom in enumerate(catoms[:-3]):
	        UB[catoms[i]][catoms[i+3]] = 2.8
	        UB[catoms[i+3]][catoms[i]] = 2.8        
	        LB[catoms[i]][catoms[i+3]] = 2.8
	        LB[catoms[i+3]][catoms[i]] = 2.8                              
    for i in range(natoms):
        for j in range(i):
            # fill LBs with sums of vdW radii and UBs with arbitrary large cutoff
            if LB[i][j] == 0:
                LB[i][j] = vdwrad[mol.getAtom(i).sym] + vdwrad[mol.getAtom(j).sym]
                LB[j][i] = vdwrad[mol.getAtom(i).sym] + vdwrad[mol.getAtom(j).sym]            
                UB[i][j] = 100
                UB[j][i] = 100              
    return LB,UB

## Triangle inequality bounds smoothing
#  
#  Copied from ref [2], pp. 252-253
#  
#  Scales O(N^3).
#  @param LB Lower bounds matrix
#  @param UB Upper bounds matrix
#  @param natoms Number of atoms in molecule
#  @return Triangularized bounds matrices
def Triangle(LB,UB,natoms):
    LL = LB
    UL = UB
    for k in range(natoms):
        for i in range(natoms-1):
            for j in range(i,natoms):
                if UL[i][j] > UL[i][k] + UL[k][j]:
                    UL[i][j] = UL[i][k] + UL[k][j]
                    UL[j][i] = UL[i][k] + UL[k][j]
                if LL[i][j] < LL[i][k] - UL[k][j]:
                    LL[i][j] = LL[i][k] - UL[k][j]
                    LL[j][i] = LL[i][k] - UL[k][j]
                else:
                    if LL[i][j] < LL[j][k] - UL[k][i]:
                        LL[i][j] = LL[j][k] - UL[k][i]
                        LL[j][i] = LL[j][k] - UL[k][i]
    return LL,UL

## Metrization to select random in-range distances
#  
#  Copied from ref [2], pp. 253-254
#  @param LB Lower bounds matrix
#  @param UB Upper bounds matrix
#  @param natoms Number of atoms in molecule
#  @param Full Full metrization (scales O(N^5), default false)
#  @param seed Random number seed (default none) 
#  @return Distance matrix
def Metrize(LB,UB,natoms,Full=False,seed=False):
    if seed:
        numpy.random.seed(seed)
    D = np.zeros((natoms,natoms))
    LB,UB = Triangle(LB,UB,natoms)
    for i in range(natoms-1):
        for j in range(i,natoms):
            if Full:
                LB,UB = Triangle(LB,UB,natoms)
            D[i][j] = np.random.uniform(LB[i][j],UB[i][j])
            D[j][i] = D[i][j]
    return D

## Get distances of each atom to CM given the distance matrix
#  
#  Copied from ref [2], pp. 309
#  @param D Distance matrix
#  @param natoms Number of atoms in molecule
#  @return Vector of CM distances, flag for successful search
def GetCMDists(D,natoms):
    D0 = np.zeros(natoms)
    try:
        status = True
        for i in range(natoms):
            for j in range(natoms):
                D0[i] += D[i][j]**2/natoms
            for j in range(natoms):
                for k in range(j,natoms):
                    D0[i] -= (D[j][k])**2/natoms**2
            D0[i] = sqrt(D0[i])
    except ValueError:
        status = False
    return D0,status

##  Get metric matrix from distance matrix and CM distances
#  
#  Copied from ref [1], pp. 306
#  @param D Distance matrix
#  @param D0 Vector of CM distances
#  @param natoms Number of atoms in molecule
#  @return Metric matrix
def GetMetricMatrix(D,D0,natoms):
    G = np.zeros((natoms,natoms))
    for i in range(natoms):
        for j in range(natoms):
            G[i][j] = (D0[i]**2 + D0[j]**2 - D[i][j]**2)/2
    return G  

##  Gets 3 largest eigenvalues and corresponding eigenvectors of metric matrix
#  @param G Metric matrix
#  @param natoms Number of atoms in molecule
#  @return Three largest eigenvalues and corresponding eigenvectors
def Get3Eigs(G,natoms):
    L = np.zeros((3,3))
    V = np.zeros((natoms,3))
    l,v = np.linalg.eigh(G)
    for i in [0,1,2]:
        L[i][i] = sqrt(l[natoms-1-i])
        V[:,i] = v[:,natoms-1-i]
    return L,V

## Computes distance error function for scipy optimization
#  
#  Copied from E3 in pp. 311 of ref. [1]
#  @param x 1D array of coordinates to be optimized
#  @param *args Other parameters (refer to scipy.optimize docs)
#  @return Objective function
def DistErr(x,*args):
    E = 0
    LB,UB,natoms = args
    for i in range(natoms-1):
        for j in range(i+1,natoms):
            ri = [x[3*i],x[3*i+1],x[3*i+2]]
            rj = [x[3*j],x[3*j+1],x[3*j+2]]
            dij = distance(ri,rj)
            uij = UB[i][j]
            lij = LB[i][j]            
            E += (dij**2/(uij**2) - 1)**2 
            E += (2*lij**2/(lij**2 + dij**2) - 1)**2
    return np.asarray(E)

## Computes gradient of distance error function for scipy optimization
#  
#  Copied from E3 in pp. 311 of ref. [1]
#  @param x 1D array of coordinates to be optimized
#  @param *args Other parameters (refer to scipy.optimize docs)
#  @return Objective function gradient
def DistErrGrad(x,*args):
    LB,UB,natoms = args
    g = np.zeros(3*natoms)
    for i in range(natoms):
        jr = range(natoms)
        jr.remove(i)
        for j in jr:
            ri = [x[3*i],x[3*i+1],x[3*i+2]]
            rj = [x[3*j],x[3*j+1],x[3*j+2]]
            dij = distance(ri,rj)
            uij = UB[i][j]        
            lij = LB[i][j]
            g[3*i] += (4*((dij/uij)**2-1)/(uij**2) - (8/lij**2)*(2*(lij**2/(lij**2+dij**2))-1)/((1+(dij/lij)**2)**2))*(x[3*i]-x[3*j]) # xi
            g[3*i+1] += (4*((dij/uij)**2-1)/(uij**2) - (8/lij**2)*(2*(lij**2/(lij**2+dij**2))-1)/((1+(dij/lij)**2)**2))*(x[3*i+1]-x[3*j+1]) # yi  
            g[3*i+2] += (4*((dij/uij)**2-1)/(uij**2) - (8/lij**2)*(2*(lij**2/(lij**2+dij**2))-1)/((1+(dij/lij)**2)**2))*(x[3*i+2]-x[3*j+2]) # zi
    return g

## Further cleans up with OB FF and saves to a new mol3D object
#  
#  Note that distance geometry tends to produce puckered aromatic rings because of the lack of explicit impropers, see Riniker et al. JCIM (2015) 55, 2562-74 for details.
#  
#  Hence, a FF optimization (with connection atoms constrained) is recommended to clean up the structure.
#  @param X Array of coordinates
#  @param mol mol3D of original molecule
#  @param ffclean Flag for OB FF cleanup (default True)
#  @param catoms List of connection atoms (default empty), used to generate FF constraints if specified
#  @return mol3D of new conformer
def SaveConf(X,mol,ffclean=True,catoms=[]):
    conf3D = mol3D()
    conf3D.copymol3D(mol)
    # set coordinates using OBMol to keep bonding info
    OBMol = conf3D.OBMol
    for i,atom in enumerate(openbabel.OBMolAtomIter(OBMol)):
        atom.SetVector(X[i,0],X[i,1],X[i,2])
    if ffclean:
        ff = openbabel.OBForceField.FindForceField('mmff94')
        constr = openbabel.OBFFConstraints()
        # constrain connecting atoms
        for catom in catoms:
            constr.AddAtomConstraint(catom+1)
        s = ff.Setup(OBMol,constr)
        if not s:
            print('FF setup failed')
        for i in range(200):
            ff.SteepestDescent(10)
            ff.ConjugateGradients(10)
        ff.GetCoordinates(OBMol)
        conf3D.OBMol = OBMol  
    conf3D.convert2mol3D()
    return conf3D

## Uses distance geometry to get a random conformer.
#  @param mol mol3D of molecule
#  @param catoms List of connection atoms (default empty), used to generate additional constraints if specified (see GetBoundsMatrices())
#  @return mol3D of new conformer
def GetConf(mol,catoms=[]):
    natoms = mol.natoms
    mol.createMolecularGraph()
    A = mol.graph
    A = A + np.dot(A,np.transpose(A))
    start = time.time()
    LB,UB = GetBoundsMatrices(mol,natoms,catoms,A)
    status = False
    while not status:
	    D = Metrize(LB,UB,natoms)
	    D0,status = GetCMDists(D,natoms)
    G = GetMetricMatrix(D,D0,natoms)
    L,V = Get3Eigs(G,natoms)
    X = np.dot(V,L) # get projection
    x = np.reshape(X,3*natoms)
    res1 = optimize.fmin_cg(DistErr,x,fprime=DistErrGrad,gtol=0.1,args=(LB,UB,natoms),disp=0)
    Opt = time.time()
    print('Optimize',str(Opt-start))
    X = np.reshape(res1,(natoms,3))
    conf3D = SaveConf(X,mol,True,catoms)
    ff = time.time()
    print('FF',str(ff-Opt))
    return conf3D

# for testing
#mol,emsg = lig_load('c1ccc(c(c1)C=NCCN=Cc2ccccc2[O-])[O-]')
#mol,emsg = lig_load('N(C)1CCN(C)CCCN(C)CCN(C)CCC1')
#catoms = [7,10,18,19]
#catoms = [0,4,9,13]
#mol.convert2mol3D()
#conf = GetConf(mol,catoms)
#conf.writexyz('conf')
