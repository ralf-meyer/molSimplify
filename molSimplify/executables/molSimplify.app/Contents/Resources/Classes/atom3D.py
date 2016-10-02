# Written by Tim Ioannidis for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
########   Defines class of 3D atoms that    #############
########     will be used to manipulate      #############
########   coordinates and other properties  #############
##########################################################

# import modules
from math import sqrt 
from Classes.globalvars import globalvars


class atom3D:
    ################################
    ### constructor for 3D atoms ###
    ################################
    def __init__(self,Sym='C',xyz=[0.0,0.0,0.0]): 
        # INPUT
        #   - Sym: element symbol
        #   - xyz: list with 3d point
        """ Create a new atom object """
        self.sym = Sym
        globs = globalvars()
        amass = globs.amass()
        if Sym not in amass: # assign default values if not in dictionary
            print("We didn't find the atomic mass of %s in the dictionary. Assigning default value of 12!\n" %(Sym))
            self.mass = 12 # default atomic mass
            self.atno = 6 # default atomic number
            self.rad = 0.75 # default atomic radius
        else:
            self.mass = amass[Sym][0] # atomic mass
            self.atno = amass[Sym][1] # atomic number
            self.rad = amass[Sym][2] # atomic covalent radius
        self.__xyz = xyz # coords
            
    ################################
    ### get coordinates of atoms ###
    ################################
    def coords(self): 
        # OUTPUT
        #   list with xyz coordinates
        x,y,z = self.__xyz
        return [x,y,z]
        
    ##########################################
    ### get distance with other atom in 3D ###
    ##########################################
    def distance(self,atom2):
        # INPUT
        #   - atom2: second atom
        # OUTPUT
        #   distance between the two atoms
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return sqrt(dx*dx+dy*dy+dz*dz)
        
    #############################################
    ### get distance vector with another atom ###
    #############################################
    def distancev(self,atom2):
        # INPUT
        #   - atom2: second atom
        # OUTPUT
        #   distance vector between the two atoms
        xyz = self.coords()
        point = atom2.coords()
        dx = xyz[0]-point[0]
        dy = xyz[1]-point[1]
        dz = xyz[2]-point[2]
        return [dx,dy,dz]
        
    ################################
    ### check if atom is a metal ###
    ################################
    def ismetal(self):
        # OUTPUT
        #   flag for metal or not
        if self.sym in globalvars().metals():
            return True
        else:
            return False
    ########################################
    ### set coordinates of atom in space ###
    ########################################
    def setcoords(self,xyz):
        # INPUT
        #   - xyz: list with xyz coordinates
        self.__xyz[0] = xyz[0]
        self.__xyz[1] = xyz[1]
        self.__xyz[2] = xyz[2]
        
    ##########################
    ### get symbol of atom ###
    ##########################
    def symbol(self):
        # OUTPUT
        #   element symbol for atom
        return self.sym
    
    #########################################
    ### move atom in space by vector dxyz ###
    #########################################
    def translate(self,dxyz):
        # INPUT
        #   - dxyz: translation vector
        x,y,z = self.__xyz
        self.__xyz[0] = x + dxyz[0]
        self.__xyz[1] = y + dxyz[1]
        self.__xyz[2] = z + dxyz[2]
        
    ######################################
    ### print methods for atom3D class ###
    ######################################
    def __repr__(self):
        # OUTPUT
        #   - ss: string with methods associated with the class
        """ when calls mol3D object without attribute e.g. t """
        ss = "\nClass atom3D has the following methods:\n"
        for method in dir(self):
            if callable(getattr(self, method)):
                ss += method +'\n'
        return ss
