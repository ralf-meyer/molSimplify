# @file AA3D.py
#  Defines AA3D class and contains useful manipulation/retrieval routines.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

# imports
import os
from math import sqrt

class AA3D:
    """Holds information about an amino acid, used to do manipulations.  Reads information from structure file (pdb, cif) or is directly built from molsimplify.
    
    """
    
    def __init__(self, three_lc='GLY', chain='undef', id=-1, occup=1.00):
        # List of atom3D objects
        self.atoms = []
        # Number of atoms
        self.natoms = 0
        # 3-letter code of amino acid (in all caps)
        self.three_lc = three_lc # if no name is specified, defaults to glycine
        # Mass of molecule
        self.mass = 0
        # Size of molecule
        self.size = 0
        # Chain of amino acid
        self.chain = chain
        # ID of amino acid (position in chain)
        self.id = id
        # Occupancy of amino acid in chain
        self.occup = occup
        # Peptide bond atoms
        self.c = None
        self.n = None
        # Bonds
        self.bonds = {}
            
    def identify(self):
        """ States whether the amino acid is (positively/negatively) charged, polar, or hydrophobic.
        
        Returns
        -------
        aa_type : string
            Positively charged, Negatively charged, Polar, Hydrophobic
                
        """
        if self.three_lc == "ARG" or self.three_lc == "LYS":
            return "Positively charged"
        elif self.three_lc == "ASP" or self.three_lc == "GLU":
            return "Negatively charged"
        elif self.three_lc in {"GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS"}:
            return "Polar"
        else:
            return "Hydrophobic"
            
    def getGreek(self, greek):
        """ Finds the Greek lettered carbon(s) or other atom(s) of the user's choice.
        
        Parameters
        ----------
        greek : string
            The Greek lettered atom (e.g. alpha carbon) we want.  Inputs should be form 'CA' or similar.
                
        Returns
        -------
        greek_atoms : list of atom3Ds
            A list of atom3D class objects that contains the Greek lettered atom(s) we want.
                
        """
        greek_atoms = []
        for a in self.atoms:
                if greek in a.greek: greek_atoms.append(a)
        return greek_atoms

    def coords(self):
        """Method to obtain string of coordinates in amino acid.

        Returns
        -------
            coord_string : string
                String of molecular coordinates with atom identities in XYZ format.
        """
        coord_string = ''  # initialize returning string
        coord_string += "%d \n\n" % self.natoms
        for atom in self.atoms:
            xyz = atom.coords()
            coord_string += "%s \t%f\t%f\t%f\n" % (atom.sym, xyz[0], xyz[1], xyz[2])
        return coord_string
            
    def centermass(self):
        """Computes coordinates of center of mass of amino acid.
        Returns
        -------
            center_of_mass : list
                Coordinates of center of mass. List of length 3: (X, Y, Z).
        """

        center_of_mass = [0, 0, 0]
        mmass = 0
        # loop over atoms in molecule
        if self.natoms > 0:
            for atom in self.atoms:
                # calculate center of mass (relative weight according to atomic mass)
                xyz = atom.coords()
                center_of_mass[0] += xyz[0] * atom.mass
                center_of_mass[1] += xyz[1] * atom.mass
                center_of_mass[2] += xyz[2] * atom.mass
                mmass += atom.mass
            # normalize
            center_of_mass[0] /= mmass
            center_of_mass[1] /= mmass
            center_of_mass[2] /= mmass
        else:
            center_of_mass = False
            print(
                'ERROR: Center of mass calculation failed. Structure will be inaccurate.\n')
        return center_of_mass

    def getBonds(self):
        """ Gets the bonds between atoms within the amino acid.

        Returns
        -------
            bonds :  dictionary (represents an undirected graph)
                Keyed by atom3D atoms in the amino acid
                Valued by a set consisting of bonded atoms
        """
        for key in self.atoms:
            if key not in self.bonds.keys():
                self.bonds[key] = set()
                if key.greek == 'C':
                    for a in self.atoms:
                        if a.greek == 'O' or a.greek == "CA":
                            self.bonds[key].add(a)
                            self.bonds[a].add(key)
                if key.greek == 'CA':
                    for a in self.atoms:
                        if a.greek == 'CB' or a.greek == "N":
                            self.bonds[key].add(a)
                            self.bonds[a].add(key)
                if "G" in key.greek:
                    for a in self.atoms:
                        if 'B' in a.greek or "D" in a.greek:
                            self.bonds[key].add(a)
                            self.bonds[a].add(key)
                if "E" in key.greek:
                    for a in self.atoms:
                        if 'Z' in a.greek or "D" in a.greek:
                            self.bonds[key].add(a)
                            self.bonds[a].add(key)
        return bonds

    def setPeptideAtoms(self):
        """ Sets the atoms involved in peptide bonding as attributes.

        """
        for a in self.atoms:
            if a.greek == "C":
                self.c = a
            elif a.greek == "N":
                self.n = a
                    
                            
            

