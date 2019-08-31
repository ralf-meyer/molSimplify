# @file atom3D.py
#  Defines atom3D class and contains useful manipulation/retrieval routines.
#
#  Written by Tim Ioannidis and JP Janet for HJK Group
#
#  Dpt of Chemical Engineering, MIT

from math import sqrt
from molSimplify.Classes.globalvars import globalvars

# Class for atoms that will be used to manipulate coordinates and other properties


class atom3D:
    # Constructor
    #  @param self The object pointer
    #  @param Sym Element symbol
    #  @param xyz List of coordinates
    #  @param name optional unique identifier
    def __init__(self, Sym='C', xyz=[0.0, 0.0, 0.0], name=False):
                # Element symbol
        self.sym = Sym
        globs = globalvars()
        amass = globs.amass()
        if Sym not in amass:  # assign default values if not in dictionary
            print(("We didn't find the atomic mass of %s in the dictionary. Assigning default value of 12!\n" % (Sym)))
            # Atomic mass
            self.mass = 12  # default atomic mass
            # Atomic number
            self.atno = 6  # default atomic number
            # Covalent radius
            self.rad = 0.75  # default atomic radius
        else:
            self.mass = amass[Sym][0]
            self.atno = amass[Sym][1]
            self.rad = amass[Sym][2]
        # Flag for freezing in optimization
        self.frozen = False
        # Flag for atom name
        if name:
            self.name = name
        else:
            self.name = Sym
        # flag for metal
        self.metal = None

        # Coordinates
        self.__xyz = xyz

    # Get coordinates
    #  @param self The object pointer
    #  @return List of coordinates
    def coords(self):
        x, y, z = self.__xyz
        return [x, y, z]

    # Get distance with another atom
    #  @param self The object pointer
    #  @param atom2 atom3D of second atom
    #  @return Distance scalar
    def distance(self, atom2):
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return sqrt(dx*dx+dy*dy+dz*dz)

    # Get distance vector with another atom
    #  @param self The object pointer
    #  @param atom2 atom3D of second atom
    #  @return Distance vector
    def distancev(self, atom2):
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return [dx, dy, dz]

    # Check if atom is a metal
    #  @param self The object pointer
    #  @return ismetal bool
    def ismetal(self):
        if self.metal is None:
            if self.sym in globalvars().metals():
                self.metal = True
            else:
                self.metal = False
        return self.metal

    # Set 3D coordinates
    #  @param self The object pointer
    #  @param xyz Coordinate vector
    def setcoords(self, xyz):
        self.__xyz[0] = xyz[0]
        self.__xyz[1] = xyz[1]
        self.__xyz[2] = xyz[2]

    # Get symbol
    #  @param self The object pointer
    #  @return Element symbol
    def symbol(self):
        return self.sym

    # mutate
    #  @param self The object pointer
    #  @parm newType the new atom type
    def mutate(self, newType):
        globs = globalvars()
        amass = globs.amass()
        if newType not in list(amass.keys()):
            print(('Error, unknown atom atom type transformation to ' + str(newType)))
            print('no changes made')
        else:
            self.mass = amass[newType][0]
            self.atno = amass[newType][1]
            self.rad = amass[newType][2]
            self.name = newType
            self.sym = newType

    # Translate atom
    #  @param self The object pointer
    #  @param dxyz Distance vector to translate
    def translate(self, dxyz):
        x, y, z = self.__xyz
        self.__xyz[0] = x + dxyz[0]
        self.__xyz[1] = y + dxyz[1]
        self.__xyz[2] = z + dxyz[2]

    # Print methods
    #  @param self The object pointer
    #  @return String with methods
    def __repr__(self):
        ss = "\nClass atom3D has the following methods:\n"
        for method in dir(self):
            if callable(getattr(self, method)):
                ss += method + '\n'
        return ss
