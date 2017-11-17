## @file geometry.py
#  Contains many useful 3D Euclidean geometric manipulation routines.
#
#  Unless otherwise stated, all "points" refer to 3-element lists.
#  
#  Written by Tim Ioannidis for HJK Group
#  
#  Dpt of Chemical Engineering, MIT

import sys
import copy
from numpy import arccos, cross, dot, pi, transpose
from numpy import sin, cos, mat, array, arctan2
from numpy.linalg import det, svd
from math import pi ,sin, cos, sqrt
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.atom3D import atom3D

## Euclidean norm
#  @param u Vector
#  @return Norm of u
def norm(u):
    d = 0.0
    for u0 in u:
        d += (u0*u0)
    d = sqrt(d)
    return d
    
## Normalize a vector
#  @param u Vector
#  @return Normalized vector
def normalize(u):
    d = norm(u)
    un = []
    if d > 1.0e-13:
        un.append(u/d)
    return un

## Euclidean distance between points
#  @param R1 Point 1
#  @param R2 Point 2
#  @return Euclidean distance
def distance(R1,R2):
    dx = R1[0] - R2[0] 
    dy = R1[1] - R2[1] 
    dz = R1[2] - R2[2] 
    d = sqrt(dx**2+dy**2+dz**2)
    return d

## Element-wise vector difference
#  @param r1 Point 1
#  @param r2 Point 2
#  @return Vector difference
def vecdiff(r1,r2):
    dr = [a-b for a,b in zip(r1,r2)]
    return dr
    
## Vector midpoint
#  @param r1 Point 1
#  @param r2 Point 2
#  @return Vector midpoint
def midpt(r1,r2):
    m = [0.5*(a+b) for a,b in zip(r1,r2)]
    return m   
    
## Checks if three points are collinear
#  @param R1 Point 1
#  @param R2 Point 2
#  @param R3 Point 3    
#  @return Collinear flag
def checkcolinear(R1,R2,R3):
    dr1 = vecdiff(R2,R1)
    dr2 = vecdiff(R1,R3)
    dd = cross(array(dr1),array(dr2))
    if norm(dd) < 1.e-01:
        return True
    else:
        return False
        
## Checks if four points are coplanar
#  @param R1 Point 1
#  @param R2 Point 2
#  @param R3 Point 3    
#  @param R4 Point 4
#  @return Coplanar flag
def checkplanar(R1,R2,R3,R4):
    r31 = vecdiff(R3,R1)
    r21 = vecdiff(R2,R1)
    r43 = vecdiff(R4,R3)
    cr0 = cross(array(r21),array(r43))
    dd = dot(r31,cr0)
    if abs(dd) < 1.e-1:
        return True
    else:
        return False

## Computes angle between two vectors
#  @param r1 Vector 1
#  @param r2 Vector 2
#  @return Angle between vectors in degrees
def vecangle(r1,r2):
    if(norm(r2)*norm(r1) > 1e-16):
        theta = 180*arccos(dot(r2,r1)/(norm(r2)*norm(r1)))/pi
    else:
        theta = 0.0
    return theta

## Gets point given reference point, direction vector and distance
#  @param Rr Reference point
#  @param dist Distance
#  @param u Direction vector
#  @return Final point
def getPointu(Rr, dist, u):
    # get float bond length
    bl = float(dist)
    # get unit vector through line r = r0 + t*u
    t = bl/norm(u) # get t as t=bl/norm(r1-r0)
    # get point
    P = [0,0,0]
    P[0] = t*u[0]+Rr[0]
    P[1] = t*u[1]+Rr[1]
    P[2] = t*u[2]+Rr[2]
    return P

## Gets angle between three points (r10 and r21) and and the normal vector to the plane containing three points
#  @param r0 Point 0
#  @param r1 Point 1
#  @param r2 Point 2
#  @return Angle in degrees
#  @return Normal vector
def rotation_params(r0,r1,r2):
    r10 = [a-b for a,b in zip(r1,r0)]
    r21 = [a-b for a,b in zip(r2,r1)]
    # angle between r10 and r21
    if(norm(r21)*norm(r10) > 1e-16):
        theta = 180*arccos(dot(r21,r10)/(norm(r21)*norm(r10)))/pi
    else:
        theta = 0.0
    # get normal vector to plane r0 r1 r2
    u = cross(r21,r10)
    # check for collinear case
    if norm(u) < 1e-16:
        # pick random perpendicular vector
        if (abs(r21[0]) > 1e-16):
            u = [(-r21[1]-r21[2])/r21[0],1,1]
        elif (abs(r21[1]) > 1e-16):
            u = [1,(-r21[0]-r21[2])/r21[1],1]
        elif (abs(r21[2]) > 1e-16):
            u = [1,1,(-r21[0]-r21[1])/r21[2]]
    return theta,u

## Aligns (translates and rotates) two molecules to minimize RMSD using the Kabsch algorithm
#  @param mol0 mol3D of molecule to be aligned
#  @param mol1 mol3D of reference molecule
#  @return Aligned mol3D
def kabsch(mol0,mol1):
    # translate to align centroids with origin
    mol0 = setPdistance(mol0,mol0.centersym(),[0,0,0],0)
    mol1 = setPdistance(mol1,mol1.centersym(),[0,0,0],0)    
    # get coordinates and matrices P,Q
    P, Q = [],[]
    for atom0,atom1 in zip(mol0.getAtoms(),mol1.getAtoms()):
        P.append(atom0.coords())
        Q.append(atom1.coords())
    # Computation of the covariance matrix
    C = dot(transpose(P), Q)
    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = svd(C)
    d = (det(V) * det(W)) < 0.0
    # Create Rotation matrix U
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = dot(V, W)
    # Rotate P
    P = dot(P, U)
    # write back coordinates
    for i,atom in enumerate(mol0.getAtoms()):
        atom.setcoords(P[i])
    return mol0
    
## Reflects point about plane defined by its normal vector and a point on the plane
#  @param u Normal vector to plane
#  @param r Point to be reflected
#  @param Rp Reference point on plane
#  @return Reflected point
def ReflectPlane(u,r,Rp):
    un = norm(u)
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    # construct augmented vector rr = [r;1]
    d = -u[0]*Rp[0]-u[1]*Rp[1]-u[2]*Rp[2]
    # reflection matrix
    R=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    rn = [0,0,0]
    R[0][0] = 1-2*u[0]*u[0]
    R[0][1] = -2*u[0]*u[1] 
    R[0][2] = -2*u[0]*u[2] 
    R[0][3] = -2*u[0]*d
    R[1][0] = -2*u[1]*u[0] 
    R[1][1] = 1-2*u[1]*u[1] 
    R[1][2] = -2*u[1]*u[2] 
    R[1][3] = -2*u[1]*d
    R[2][0] = -2*u[2]*u[0]
    R[2][1] = -2*u[1]*u[2]
    R[2][2] = 1-2*u[2]*u[2] 
    R[2][3] = -2*u[2]*d
    R[3][3] = 1 
    # get new point
    rn[0] = R[0][0]*r[0]+R[0][1]*r[1]+R[0][2]*r[2] + R[0][3] 
    rn[1] = R[1][0]*r[0]+R[1][1]*r[1]+R[1][2]*r[2] + R[1][3] 
    rn[2] = R[2][0]*r[0]+R[2][1]*r[1]+R[2][2]*r[2] + R[2][3] 
    return rn

## Rotates point about axis defined by direction vector and point on axis
#  @param u Direction vector of axis
#  @param rp Reference point along axis
#  @param r Point to be rotated
#  @param theta Angle of rotation in RADIANS
#  @return Rotated point
def PointRotateAxis(u,rp,r,theta):
    # construct augmented vector rr = [r;1]
    rr = r
    rr.append(1)
    # rotation matrix about arbitrary line through rp
    R=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    rn = [0,0,0]
    R[0][0] = cos(theta)+u[0]**2*(1-cos(theta))
    R[0][1] = u[0]*u[1]*(1-cos(theta))-u[2]*sin(theta)
    R[0][2] = u[0]*u[2]*(1-cos(theta))+u[1]*sin(theta)
    R[0][3] = (rp[0]*(u[1]**2+u[2]**2)-u[0]*(rp[1]*u[1]+rp[2]*u[2]))*(1-cos(theta))
    R[0][3] += (rp[1]*u[2]-rp[2]*u[1])*sin(theta)
    R[1][0] = u[1]*u[0]*(1-cos(theta))+u[2]*sin(theta)
    R[1][1] = cos(theta)+u[1]**2*(1-cos(theta))
    R[1][2] = u[1]*u[2]*(1-cos(theta))-u[0]*sin(theta)
    R[1][3] = (rp[1]*(u[0]**2+u[2]**2)-u[1]*(rp[0]*u[0]+rp[2]*u[2]))*(1-cos(theta))
    R[1][3] += (rp[2]*u[0]-rp[0]*u[2])*sin(theta)
    R[2][0] = u[2]*u[0]*(1-cos(theta))-u[1]*sin(theta)
    R[2][1] = u[2]*u[1]*(1-cos(theta))+u[0]*sin(theta)
    R[2][2] = cos(theta)+u[2]**2*(1-cos(theta))
    R[2][3] = (rp[2]*(u[0]**2+u[1]**2)-u[2]*(rp[0]*u[0]+rp[1]*u[1]))*(1-cos(theta))
    R[2][3] += (rp[0]*u[1]-rp[1]*u[0])*sin(theta)
    R[3][3] = 1
    # get new point
    rn[0] = R[0][0]*r[0]+R[0][1]*r[1]+R[0][2]*r[2] + R[0][3]
    rn[1] = R[1][0]*r[0]+R[1][1]*r[1]+R[1][2]*r[2] + R[1][3]
    rn[2] = R[2][0]*r[0]+R[2][1]*r[1]+R[2][2]*r[2] + R[2][3]
    return rn
    
## Translates point in spherical coordinates
#  @param Rp Origin of sphere
#  @param p0 Point to be translated
#  @param D [final radial distance, change in polar phi, change in azimuthal theta] in RADIANS
#  @return Translated point
def PointTranslateSph(Rp,p0,D):    
    # translate to origin
    ps=[0,0,0]
    ps[0] = p0[0] - Rp[0]
    ps[1] = p0[1] - Rp[1]
    ps[2] = p0[2] - Rp[2]  
    # get initial spherical coords
    r0 = norm(ps) 
    if (r0 < 1e-16):
        phi0 = 0.5*pi
        theta0 = 0 
    else :
        phi0 = arccos(ps[2]/r0) # z/r
        theta0 = arctan2(ps[1],ps[0]) # y/x
    # get new point
    p = [0,0,0]
    p[0] = (D[0])*sin(phi0+D[2])*cos(theta0+D[1]) + Rp[0]
    p[1] = (D[0])*sin(phi0+D[2])*sin(theta0+D[1]) + Rp[1]
    p[2] = (D[0])*cos(phi0+D[2]) + Rp[2]
    return p

## Converts spherical translation vector into Cartesian translation vector
#  @param Rp Origin of sphere
#  @param p0 Point to be translated
#  @param D [final radial distance, change in polar phi, change in azimuthal theta] in RADIANS
#  @return Translation vector
def PointTranslatetoPSph(Rp,p0,D):    
    # translate to origin
    ps=[0,0,0]
    ps[0] = p0[0] - Rp[0]
    ps[1] = p0[1] - Rp[1]
    ps[2] = p0[2] - Rp[2]  
    # get current spherical coords
    r0 = norm(ps) 
    if (r0 < 1e-16):
        phi0 = 0.5*pi
        theta0 = 0 
    else :
        phi0 = arccos(ps[2]/r0) # z/r
        theta0 = arctan2(ps[1],ps[0]) # y/x
    # get translation vector
    p = [0,0,0]
    p[0] = D[0]*sin(phi0+D[2])*cos(theta0+D[1]) 
    p[1] = D[0]*sin(phi0+D[2])*sin(theta0+D[1]) 
    p[2] = D[0]*cos(phi0+D[2]) 
    return p

## Rotates point about Cartesian axes defined relative to given origin
#  @param Rp Cartesian origin
#  @param p0 Point to be rotated
#  @param D [theta-x, theta-y, theta-z] in RADIANS
#  @return Rotated point
def PointRotateSph(Rp,p0,D):
    # translate to origin (reference)
    ps=[0,0,0]
    ps[0] = p0[0] - Rp[0]
    ps[1] = p0[1] - Rp[1]
    ps[2] = p0[2] - Rp[2]  
    # build 3D rotation matrices about x,y,z axes
    Mx=[[1, 0, 0],[0, cos(D[0]), -sin(D[0])],[0, sin(D[0]), cos(D[0])]]
    My=[[cos(D[1]), 0, sin(D[1])],[0, 1, 0],[-sin(D[1]), 0, cos(D[1])]]
    Mz=[[cos(D[2]), -sin(D[2]), 0],[sin(D[2]), cos(D[2]), 0],[0, 0, 1]]
    # get full rotation matrix
    M = array(mat(Mx)*mat(My)*mat(Mz))
    p=[0.0, 0.0, 0.0]
    # rotate atom and translate it back from origin
    p[0] = M[0][0]*ps[0] + M[0][1]*ps[1] + M[0][2]*ps[2] + Rp[0]
    p[1] = M[1][0]*ps[0] + M[1][1]*ps[1] + M[1][2]*ps[2] + Rp[1]
    p[2] = M[2][0]*ps[0] + M[2][1]*ps[1] + M[2][2]*ps[2] + Rp[2]
    return p
    
## Reflects molecule about plane defined by its normal vector and a point on the plane
#
#  Loops over ReflectPlane().
#  @param mol mol3D of molecule to be reflected
#  @param u Normal vector to plane
#  @param Rp Reference point on plane
#  @return mol3D of reflected molecule
def reflect_through_plane(mol,u,Rp):
    un = norm(u)
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = ReflectPlane(u,atom.coords(),Rp)
        atom.setcoords(Rt)
    return mol

## Rotates molecule about axis defined by direction vector and point on axis
#
#  Loops over PointRotateAxis().
#  @param mol mol3D of molecule to be rotated
#  @param Rp Reference point along axis
#  @param u Direction vector of axis
#  @param theta Angle of rotation in DEGREES
#  @return mol3D of rotated molecule
def rotate_around_axis(mol,Rp,u,theta):
    un = norm(u)
    theta = (theta/180.0)*pi
    if (un > 1e-16):
        u[0] = u[0]/un
        u[1] = u[1]/un
        u[2] = u[2]/un
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateAxis(u,Rp,atom.coords(),theta)
        atom.setcoords(Rt)
    return mol

## Translates molecule such that a given point in the molecule is at a given distance from a reference point
#  
#  The molecule is moved along the axis given by the two points.
#  @param mol mol3D of molecule to be translated
#  @param Rr Point in molecule to be aligned
#  @param Rp Reference alignment point
#  @param bond Final distance of aligned point to alignment point
#  @return mol3D of translated molecule
def setPdistance(mol, Rr, Rp, bond):
    # get float bond length
    bl = float(bond)
    # get center of mass
    # get unit vector through line r = r0 + t*u
    u = [a-b for a,b in zip(Rr,Rp)]
    t = bl/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for centermass
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-Rr[0]
    dxyz[1] = Rp[1]+t*u[1]-Rr[1]
    dxyz[2] = Rp[2]+t*u[2]-Rr[2]
    # translate molecule
    mol.translate(dxyz)
    return mol
    
## Translates molecule such that a given point in the molecule is at a given distance from a reference point
#  
#  The molecule is moved along an arbitrary axis.
#  @param mol mol3D of molecule to be translated
#  @param Rr Point in molecule to be aligned
#  @param Rp Reference alignment point
#  @param bond Final distance of aligned point to alignment point
#  @param u Direction vector of axis
#  @return mol3D of translated molecule
def setPdistanceu(mol, Rr, Rp, bond, u):
    # get float bond length
    bl = float(bond)
    # get unit vector through line r = r0 + t*u
    t = bl/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for centermass
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-Rr[0]
    dxyz[1] = Rp[1]+t*u[1]-Rr[1]
    dxyz[2] = Rp[2]+t*u[2]-Rr[2]
    # translate molecule
    mol.translate(dxyz)
    return mol
    
## Translates molecule such that its center of mass is at a given distance from a reference point
#  
#  The molecule is moved along the axis given by the two points.
#  @param mol mol3D of molecule to be translated
#  @param Rp Reference alignment point
#  @param bond Final distance of aligned point to alignment point
#  @return mol3D of translated molecule
def setcmdistance(mol, Rp, bond):
    # get float bond length
    bl = float(bond)
    # get center of mass
    cm = mol.centermass()
    # get unit vector through line r = r0 + t*u
    u = [a-b for a,b in zip(cm,Rp)]
    t = bl/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for centermass
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-cm[0]
    dxyz[1] = Rp[1]+t*u[1]-cm[1]
    dxyz[2] = Rp[2]+t*u[2]-cm[2]
    # translate molecule
    mol.translate(dxyz)
    return mol

## Translates molecule in spherical coordinates based on center of mass reference
#
#  Loops over PointTranslateSph().
#  @param mol mol3D of molecule to be translated
#  @param Rr Origin of sphere
#  @param D [final radial distance, change in polar phi, change in azimuthal theta] in RADIANS
#  @return mol3D of translated molecule
def protate(mol, Rr, D):
    # convert to rad
    D[0] = float(D[0])
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # rotate/translate about reference point
    # get center of mass
    pmc = mol.centermass()
    # get translation vector that corresponds to new coords
    Rt = PointTranslateSph(Rr,pmc,D)
    # translate molecule
    mol.translate(Rt)
    return mol


## Translates molecule in spherical coordinates based on arbitrary reference
#
#  Loops over PointTranslateSph().
#  @param mol mol3D of molecule to be translated
#  @param Rr Origin of sphere
#  @param Rref Reference point in molecule
#  @param D [final radial distance, change in polar phi, change in azimuthal theta] in RADIANS
#  @return mol3D of translated molecule
def protateref(mol, Rr, Rref, D):
    # rotate/translate about reference point
    # convert to rad
    D[0] = float(D[0])
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # rotate/translate about reference point
    # get translation vector that corresponds to new coords
    Rt = PointTranslateSph(Rr,Rref,D)
    # translate molecule
    mol.translate(Rt)
    return mol

## Rotates molecule about its center of mass
#
#  Loops over PointRotateSph().
#  @param mol mol3D of molecule to be rotated
#  @param D [theta-x, theta-y, theta-z] in RADIANS
#  @return mol3D of rotated molecule
def cmrotate(mol, D):
    # convert to rad
    D[0] = (float(D[0])/180.0)*pi
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # perform rotation
    pmc = mol.centermass()
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(pmc,atom.coords(),D)
        atom.setcoords(Rt)
    return mol

## Rotates molecule about an arbitrary point
#
#  Loops over PointRotateSph().
#  @param mol mol3D of molecule to be rotated
#  @param Ref Reference point
#  @param D [theta-x, theta-y, theta-z] in RADIANS
#  @return mol3D of rotated molecule
def rotateRef(mol, Ref, D):
    # convert to rad
    D[0] = (float(D[0])/180.0)*pi
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # perform rotation
    pmc = mol.centermass()
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(Ref,atom.coords(),D)
        atom.setcoords(Rt)
    return mol
    
## Translates molecule to align point to axis at constant distance
#  @param mol mol3D of molecule to be translated
#  @param Rr Point to be aligned
#  @param Rp Reference point on axis
#  @param u Direction vector of axis
#  @return mol3D of translated molecule
def aligntoaxis(mol,Rr,Rp,u):
    # INPUT
    #   - Rr: point to be aligned
    #   - mol: molecule to be manipulated
    #   - Rp: reference point through axis
    #   - u: target axis for alignment
    # OUTPUT
    #   - mol: aligned molecule
    # get current distance
    d0 = distance(Rp,Rr)
    # normalize u
    t =d0/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for point
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-Rr[0]
    dxyz[1] = Rp[1]+t*u[1]-Rr[1]
    dxyz[2] = Rp[2]+t*u[2]-Rr[2]
    # translate molecule
    mol.translate(dxyz)
    return mol

## Translates molecule to align point to axis at arbitrary distance
#  @param mol mol3D of molecule to be translated
#  @param Rr Point to be aligned
#  @param Rp Reference point on axis
#  @param u Direction vector of axis
#  @param d Final distance from aligned point to axis
#  @return mol3D of translated molecule
def aligntoaxis2(mol,Rr,Rp,u,d):
    # normalize u
    t =d/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for point
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]-Rr[0]
    dxyz[1] = Rp[1]+t*u[1]-Rr[1]
    dxyz[2] = Rp[2]+t*u[2]-Rr[2]
    # translate molecule
    mol.translate(dxyz)
    return mol
    
## Translates point and aligns to axis
#  @param Rr Point to be aligned
#  @param Rp Reference point on axis
#  @param u Direction vector of axis
#  @param d Final distance from aligned point to axis
#  @return Translation vector
def alignPtoaxis(Rr,Rp,u,d):
    # INPUT
    #   - Rr: point to be aligned
    #   - Rp: reference point through axis
    #   - u: target axis for alignment
    #   - d: final distance target
    # OUTPUT
    #   - dxyz: final coordinates
    # normalize u
    t =d/norm(u) # get t as t=bl/norm(r1-r0)
    # get shift for point
    dxyz = [0,0,0]
    dxyz[0] = Rp[0]+t*u[0]
    dxyz[1] = Rp[1]+t*u[1]
    dxyz[2] = Rp[2]+t*u[2]
    return dxyz
        
## Rotates molecule about Cartesian axes defined relative to given origin
#
#  Loops over PointRotateSph().
#  @param mol mol3D of molecule to be rotated
#  @param Rp Cartesian origin
#  @param D [theta-x, theta-y, theta-z] in DEGREES
#  @return mol3D of rotated molecule
def pmrotate(mol, Rp, D):
    # convert to rad
    D[0] = (float(D[0])/180.0)*pi
    D[1] = (float(D[1])/180.0)*pi
    D[2] = (float(D[2])/180.0)*pi
    # perform rotation
    for atom in mol.atoms:
        # Get new point after rotation
        Rt = PointRotateSph(Rp,atom.coords(),D)
        atom.setcoords(Rt)
    return mol   
