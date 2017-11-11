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
#  @param mol mol3D of molecule
#  @param natoms Number of atoms in molecule
#  @return Lower and upper bounds matrices
def GetBoundsMatrices(mol,natoms):
    # We use only 1-2 (bond lengths) and 1-3 (bond angles) constraints, read from the FF-optimized initial conformer.
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
    # 
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
#  @return Distance matrix
def Metrize(LB,UB,natoms,Full=True):

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
#  @return Vector of CM distances
def GetCMDists(D,natoms):
    D0 = np.zeros(natoms)
    for i in range(natoms):
        for j in range(natoms):
            D0[i] += D[i][j]**2/natoms
        for j in range(natoms):
            for k in range(j,natoms):
                D0[i] -= (D[j][k])**2/natoms**2
        D0[i] = sqrt(D0[i])
    return D0

#  Get metric matrix from distance matrix and CM distances
#  
#  Copied from ref [1], pp. 306
#  @param D Distance matrix
#  @param D0 Vector of CM distances
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
            E += (distance(ri,rj)**2/(UB[i][j]**2) - 1)**2 
            E += (2*LB[i][j]**2/(LB[i][j]**2 + distance(ri,rj)**2) - 1)**2
    return E

## Further cleans up with OB FF and saves to a new mol3D object
#  @param X Array of coordinates
#  @param mol mol3D of original molecule
#  @param ffclean Flag for OB FF cleanup (default True)
#  @return mol3D of new conformer
def SaveConf(X,mol,ffclean=True):
    mol3Dnew = mol3D()
    mol3Dnew.copymol3D(mol)
    for i in range(mol3Dnew.natoms):
        mol3Dnew.getAtom(i).setcoords(X[i,:])
    mol3Dnew.convert2OBMol()
    if ffclean:
        OBMol = mol3Dnew.OBMol
        ff = openbabel.OBForceField.FindForceField('mmff94')
        s = ff.Setup(OBMol)
        if not s:
            print('FF setup failed')
        for i in range(10):
            ff.SteepestDescent(10)
            ff.ConjugateGradients(10)
        ff.GetCoordinates(OBMol)
        mol3Dnew.OBMol = OBMol  
        mol3Dnew.convert2mol3D()  
    return mol3Dnew

## Uses distance geometry to get a random conformer.
#  @param mol mol3D of molecule
#  @return mol3D of new conformer
def GetConf(mol):
    natoms = mol.natoms
    start = time.time()
    LB,UB = GetBoundsMatrices(mol,natoms)
    BM = time.time()
    print('Bounds ',str(BM-start))
    D = Metrize(LB,UB,natoms,Full=True)
    Met = time.time()
    print('Metrize ',str(Met-BM))
    D0 = GetCMDists(D,natoms)
    G = GetMetricMatrix(D,D0,natoms)
    L,V = Get3Eigs(G,natoms)
    X = np.dot(V,L) # get projection
    x = np.reshape(X,3*natoms)
    res1 = optimize.fmin_cg(DistErr,x,gtol=0.1,args=(LB,UB,natoms))
    Opt = time.time()
    print('Optimize ',str(Opt-Met))
    X = np.reshape(res1,(natoms,3))
    conf3D = SaveConf(X,mol,ffclean=True)
    ff = time.time()
    print('FF '+str(ff-Opt))
    return conf3D
    
mol = mol3D()
mol.getOBMol('NCCCN','smistring',ffclean=True)
mol.convert2mol3D()
conf3D = GetConf(mol)
print distance(conf3D.getAtom(0).coords(),conf3D.getAtom(4).coords())
