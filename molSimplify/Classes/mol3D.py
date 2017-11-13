## @file mol3D.py
#  Defines mol3D class and contains useful manipulation/retrieval routines.
#
#  Written by Tim Ioannidis and JP Janet for HJK Group
#
#  Dpt of Chemical Engineering, MIT

from math import sqrt
import numpy as np
from molSimplify.Classes.atom3D import atom3D
from molSimplify.Classes.globalvars import globalvars
import openbabel
import sys, time, os, subprocess, random, shutil, unicodedata, inspect, tempfile
from pkg_resources import resource_filename, Requirement
import xml.etree.ElementTree as ET
try:
    import PyQt5
    from molSimplify.Classes.miniGUI import *
    ## PyQt5 flag
    qtflag = True
except ImportError:
    qtflag = False
    pass

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

## Wrapper for executing bash commands
#  @param cmd Command to be executed
#  @return Stdout 
def mybash(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        stdout.append(line)
        if line == '' and p.poll() != None:
            break
    return ''.join(stdout)

## Class for molecules that will be used to manipulate coordinates and other properties
class mol3D:
    ## Constructor
    #  @param self The object pointer
    def __init__(self):
		## List of atom3D objects
        self.atoms = []
        ## Number of atoms
        self.natoms = 0
        ## Mass of molecule
        self.mass = 0
        ## Size of molecule
        self.size = 0
        ## Charge of molecule
        self.charge = 0
        ## Force field optimization settings
        self.ffopt = 'BA'
        ## Name of molecule
        self.name = ''
        ## Holder for openbabel molecule
        self.OBMol = False 
        ## List of connection atoms
        self.cat = [] 
        ## Denticity
        self.denticity = 0 
        ## Identifier
        self.ident = '' 
        ## Holder for global variables
        self.globs = globalvars()
        ## Holder for molecular graph
        self.graph = []

    ## Add atom to molecule
    #
    #  Added atom is appended to the end of the list.
    #  @param self The object pointer
    #  @param atom atom3D of atom to be added
    def addAtom(self,atom):
        self.atoms.append(atom)
        if atom.frozen:
            self.atoms[-1].frozen = True
        self.natoms += 1
        self.mass += atom.mass
        self.size = self.molsize()
        self.graph = []

    ## Aligns two molecules such that the coordinates of two atoms overlap.
    #
    #  Second molecules is translated relative to the first.
    #  No rotations are performed here. Use other functions for rotations.
    #  @param self The object pointer
    #  @param atom1 atom3D of reference atom in first molecule (not translated)
    #  @param atom2 atom3D of reference atom in second molecule (translated)
    def alignmol(self,atom1,atom2):
        ## get distance vector between atoms 1,2
        dv = atom2.distancev(atom1)
        self.translate(dv)

    ## Performs bond centric manipulation (same as Avogadro, stretching/squeezing bonds)
    #
    #  A submolecule is translated along the bond axis connecting it to an anchor atom.
    #
    #  Illustration: H3A-BH3 -> H3A----BH3 where B = idx1 and A = idx2
    #  @param self The object pointer
    #  @param idx1 Index of bonded atom containing submolecule to be moved
    #  @param idx2 Index of anchor atom
    #  @param d New bond length in Angstroms
    def BCM(self,idx1,idx2,d):
        bondv = self.getAtom(idx1).distancev(self.getAtom(idx2)) # 1 - 2
        # compute current bond length
        u = 0.0
        for u0 in bondv:
            u += (u0*u0)
        u = sqrt(u)
        dl = d - u # dl > 0: stretch, dl < 0: shrink
        dR = [i*(d/u - 1) for i in bondv]
        for i in self.getBondedAtoms(idx1):
            if i != idx2:
                self.getAtom(i).translate(dR)
        self.getAtom(idx1).translate(dR)

    ## Computes coordinates of center of mass of molecule
    #  @param self The object pointer
    #  @return List of center of mass coordinates
    def centermass(self):
        # OUTPUT
        #   - pcm: vector representing center of mass
        # initialize center of mass and mol mass
        pmc = [0, 0, 0]
        mmass = 0
        # loop over atoms in molecule
        if self.natoms > 0:
            for atom in self.atoms:
            # calculate center of mass (relative weight according to atomic mass)
                xyz = atom.coords()
                pmc[0] +=  xyz[0]*atom.mass
                pmc[1] +=  xyz[1]*atom.mass
                pmc[2] +=  xyz[2]*atom.mass
                mmass += atom.mass
            # normalize
            pmc[0] /= mmass
            pmc[1] /= mmass
            pmc[2] /= mmass
        else:
            pmc = False
            print 'ERROR: Center of mass calculation failed. Structure will be inaccurate.\n'
        return pmc

    ## Computes coordinates of center of symmetry of molecule
    #
    #  Identical to centermass, but not weighted by atomic masses.
    #  @param self The object pointer
    #  @return List of center of symmetry coordinates
    def centersym(self):
        # initialize center of mass and mol mass
        pmc = [0, 0, 0]
        # loop over atoms in molecule
        for atom in self.atoms:
        # calculate center of symmetry
            xyz = atom.coords()
            pmc[0] +=  xyz[0]
            pmc[1] +=  xyz[1]
            pmc[2] +=  xyz[2]
        # normalize
        pmc[0] /= self.natoms
        pmc[1] /= self.natoms
        pmc[2] /= self.natoms
        return pmc

    ## Converts OBMol to mol3D
    #
    #  Generally used after openbabel operations, such as when initializing a molecule from a file or FF optimizing it.
    #  @param self The object pointer
    def convert2mol3D(self):
        # initialize again
        self.initialize()
        # get elements dictionary
        elem = globalvars().elementsbynum()
        # loop over atoms
        for atom in openbabel.OBMolAtomIter(self.OBMol):
            # get coordinates
            pos = [atom.GetX(),atom.GetY(),atom.GetZ()]
            # get atomic symbol
            sym = elem[atom.GetAtomicNum() -1]
            # add atom to molecule
            self.addAtom(atom3D(sym,[pos[0],pos[1],pos[2]]))
            
    ## Converts mol3D to OBMol
    #
    #  Required for performing openbabel operations on a molecule, such as FF optimizations.
    #  @param self The object pointer
    def convert2OBMol(self):
        # write temp xyz
        self.writexyz('tempr.xyz')

        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("xyz")

        OBMol = openbabel.OBMol()
        obConversion.ReadFile(OBMol,'tempr.xyz')        
        
        self.OBMol = OBMol
        os.remove('tempr.xyz')
        
    ## Combines two molecules
    #
    #  Each atom in the second molecule is appended to the first while preserving orders.
    #  @param self The object pointer
    #  @param mol mol3D containing molecule to be added
    #  @return mol3D contaning combined molecule
    def combine(self,mol):
        cmol = self
        for atom in mol.atoms:
            cmol.addAtom(atom)
        cmol.graph = []
        return cmol

    ## Prints coordinates of all atoms in molecule
    #  @param self The object pointer
    #  @return String containing coordinates
    def coords(self):
        # OUTPUT
        #   - atom: string with xyz-style coordinates
        ss = '' # initialize returning string
        ss += "%d \n\n" % self.natoms
        for atom in self.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        return ss

    ## Returns coordinates of all atoms in molecule as a list of lists
    #  @param self The object pointer
    #  @return List of all atoms in molecule
    def coordsvect(self):
        ss = []
        for atom in self.atoms:
            xyz = atom.coords()
            ss.append(xyz)
        return ss

    ## Copies properties and atoms of another existing mol3D object into current mol3D object.
    #
    #  WARNING: NEVER EVER USE mol3D = mol0 to do this. It doesn't work.
    #  
    #  WARNING: ONLY USE ON A FRESH INSTANCE OF MOL3D.
    #  @param self The object pointer
    #  @param mol0 mol3D of molecule to be copied
    def copymol3D(self,mol0):
		# copy atoms
        for i,atom0 in enumerate(mol0.atoms):
            self.addAtom(atom3D(atom0.sym,atom0.coords()))
            if atom0.frozen:
                self.getAtom(i).frozen = True
        # copy other attributes
        self.cat = mol0.cat
        self.charge = mol0.charge
        self.denticity = mol0.denticity
        self.ident = mol0.ident
        self.ffopt = mol0.ffopt
        self.OBMol = mol0.OBMol
        
    ## Create molecular graph (connectivity matrix) from mol3D info
    #  @param self The object pointer
    def createMolecularGraph(self):
        index_set = range(0,self.natoms)
        A  = np.matrix(np.zeros((self.natoms,self.natoms)))
        for i in index_set:
            this_bonded_atoms = self.getBondedAtomsOct(i,debug=False)
            for j in index_set:
                if j in this_bonded_atoms:
                    A[i,j] = 1
        self.graph = A
        
    ## Deletes specific atom from molecule
    #
    #  Also updates mass and number of atoms, and resets the molecular graph.
    #  @param self The object pointer
    #  @param atomIdx Index of atom to be deleted
    def deleteatom(self,atomIdx):
        self.mass -= self.getAtom(atomIdx).mass
        self.natoms -= 1
        self.graph = []
        del(self.atoms[atomIdx])

    ## Freezes specific atom in molecule
    #
    #  This is for the FF optimization settings.
    #  @param self The object pointer
    #  @param atomIdx Index of atom to be frozen
    def freezeatom(self,atomIdx):
        # INPUT
        #   - atomIdx: index of atom to be frozen
        self.atoms[atomIdx].frozen = True
        
    ## Deletes list of atoms from molecule
    #
    #  Loops over deleteatom, starting from the largest index so ordering is preserved.
    #  @param self The object pointer
    #  @param Alist List of atom indices to be deleted
    def deleteatoms(self,Alist):
        for h in sorted(Alist,reverse=True):
            self.deleteatom(h)
            
    ## Freezes list of atoms in molecule
    #
    #  Loops over freezeatom(), starting from the largest index so ordering is preserved.
    #  @param self The object pointer
    #  @param Alist List of atom indices to be frozen
    def freezeatoms(self,Alist):
        # INPUT
        #   - Alist: list of atoms to be frozen
        for h in sorted(Alist,reverse=True):
            self.freezeatom(h)
            
    ## Deletes all hydrogens from molecule.
    #
    #  Calls deleteatoms, so ordering of heavy atoms is preserved.
    #  @param self The object pointer
    def deleteHs(self):
        hlist = []
        for i in range(self.natoms):
            if self.getAtom(i).sym=='H':
                hlist.append(i)
        self.deleteatoms(hlist)

    ## Gets distance between centers of mass of two molecules 
    #  @param self The object pointer
    #  @param mol mol3D of second molecule
    #  @return Center of mass distance
    def distance(self,mol):
        # INPUT
        #   - mol: second molecule
        # OUTPUT
        #   - pcm: distance between centers of mass
        cm0 = self.centermass()
        cm1 = mol.centermass()
        pmc = distance(cm0,cm1)
        return pmc

    ## Creates and saves an svg file of the molecule
    #
    #  Also renders it in a fake gui window if PyQt5 is installed.
    #  Copied from mGUI function.
    #  @param self The object pointer    
    #  @param filename Name of svg file
    def draw_svg(self,filename):
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("svg")
        obConversion.AddOption("i", obConversion.OUTOPTIONS, "") 
        ### return the svg with atom labels as a string
        svgstr = obConversion.WriteString(self.OBMol)
        namespace = "http://www.w3.org/2000/svg"
        ET.register_namespace("", namespace)
        tree = ET.fromstring(svgstr)
        svg = tree.find("{{{ns}}}g/{{{ns}}}svg".format(ns=namespace))
        newsvg = ET.tostring(svg).decode("utf-8")
        ### write unpacked svg file   
        fname = filename+'.svg'
        with open(fname, "w") as svg_file:
            svg_file.write(newsvg)
        if qtflag:
            ### Initialize fake gui   
            fakegui = miniGUI(sys.argv)
            ### Add the svg to the window
            fakegui.addsvg(fname)
            ### Show window
            fakegui.show()
        else:
            print('No PyQt5 found. SVG file written to directory.')

    ## Finds closest metal atom to a given atom
    #  @param self The object pointer
    #  @param atom0 Index of reference atom
    #  @return Index of closest metal atom
    def findcloseMetal(self,atom0):
        mm = False
        mindist = 1000
        for i,atom in enumerate(self.atoms):
            if atom.ismetal():
                if distance(atom.coords(),atom0.coords()) < mindist:
                    mindist = distance(atom.coords(),atom0.coords())
                    mm = i
        # if no metal, find heaviest atom
        if not mm:
            maxaw = 0
            for i,atom in enumerate(self.atoms):
                if atom.atno > maxaw:
                    mm = i
        return mm

    ## Finds metal atoms in molecule
    #  @param self The object pointer
    #  @return List of indices of metal atoms
    def findMetal(self):
        mm = []
        for i,atom in enumerate(self.atoms):
            if atom.ismetal():
                mm.append(i)
        return mm

    ## Finds atoms in molecule with given symbol
    #  @param self The object pointer
    #  @param sym Desired element symbol
    #  @return List of indices of atoms with given symbol
    def findAtomsbySymbol(self,sym):
        mm = []
        for i,atom in enumerate(self.atoms):
            if atom.sym==sym:
                mm.append(i)
        return mm

    ## Finds a submolecule within the molecule given the starting atom and the separating atom
    #
    #  Illustration: H2A-B-C-DH2 will return C-DH2 if C is the starting atom and B is the separating atom.
    #  
    #  Alternatively, if C is the starting atom and D is the separating atom, returns H2A-B-C.
    #  @param self The object pointer
    #  @param atom0 Index of starting atom
    #  @param atomN Index of separating atom
    #  @return List of indices of atoms in submolecule
    def findsubMol(self,atom0,atomN):
        subm = []
        conatoms = [atom0]
        conatoms += self.getBondedAtoms(atom0) # connected atoms to atom0
        if atomN in conatoms:
            conatoms.remove(atomN)  # check for atomN and remove
        subm += conatoms # add to submolecule
        while len(conatoms) > 0: # while list of atoms to check loop
            for atidx in subm: # loop over initial connected atoms
                if atidx != atomN: # check for separation atom
                    newcon = self.getBondedAtoms(atidx)
                    if atomN in newcon:
                        newcon.remove(atomN)
                    for newat in newcon:
                        if newat not in conatoms and newat not in subm:
                            conatoms.append(newat)
                            subm.append(newat)
                if atidx in conatoms:
                    conatoms.remove(atidx) # remove from list to check
        subm.sort()
        return subm
  
    ## Gets an atom with specified index
    #  @param self The object pointer
    #  @param idx Index of desired atom
    #  @return atom3D of desired atom
    def getAtom(self,idx):
        return self.atoms[idx]

    ## Gets atoms in molecule
    #  @param self The object pointer
    #  @return List of atoms in molecule
    def getAtoms(self):
        return self.atoms
        
    ## Gets number of unique elements in molecule
    #  @param self The object pointer
    #  @return List of symbols of unique elements in molecule
    def getAtomTypes(self):
        unique_atoms_list = list()
        for atoms in self.getAtoms():
            if atoms.symbol() not in unique_atoms_list:
                        unique_atoms_list.append(atoms.symbol())
        return unique_atoms_list
        
        
    ## Gets coordinates of atom with specified index
    #  @param self The object pointer
    #  @param idx Index of desired atom
    #  @return List of coordinates of desired atom
    def getAtomCoords(self,idx):
        return self.atoms[idx].coords()

    ## Gets atoms bonded to a specific atom
    #
    #  This is determined based on element-specific distance cutoffs, rather than predefined valences.
    #  
    #  This method is ideal for metals because bond orders are ill-defined.
    #
    #  For pure organics, the OBMol class provides better functionality.
    #  @param self The object pointer
    #  @param ind Index of reference atom
    #  @return List of indices of bonded atoms
    def getBondedAtoms(self,ind):
        ratom = self.getAtom(ind)
        # calculates adjacent number of atoms
        nats = []
        for i,atom in enumerate(self.atoms):
            d = distance(ratom.coords(),atom.coords())
            distance_max = 1.30*(atom.rad+ratom.rad) 
            if atom.symbol() == "C" and not ratom.symbol() == "H":
                distance_max = min(2.8,distance_max)
            if ratom.symbol() == "C" and not atom.symbol() == "H":
                distance_max = min(2.8,distance_max)
            if ratom.symbol() == "H" and atom.ismetal:
                ## tight cutoff for metal-H bonds
                distance_max = 1.1*(atom.rad+ratom.rad) 
            if atom.symbol() == "H" and ratom.ismetal:
                ## tight cutoff for metal-H bonds
                distance_max = 1.1*(atom.rad+ratom.rad) 
            if (d < distance_max and i!=ind):
                    nats.append(i)
        return nats
        
        
    ## Gets atoms bonded to a specific atom specialized for octahedral complexes
    #
    #  More sophisticated version of getBondedAtoms(), written by JP.
    #  
    #  This method specifically forbids "intruder" C and H atoms that would otherwise be within the distance cutoff in tightly bound complexes.
    #
    #  It also limits bonding atoms to the CN closest atoms (CN = coordination number).
    #  @param self The object pointer
    #  @param ind Index of reference atom
    #  @param CN Coordination number of reference atom (default 6)
    #  @param debug Debug flag (default False)
    #  @return List of indices of bonded atoms
    def getBondedAtomsOct(self,ind,CN=6,debug=False):
        # INPUT
        #   - ind: index of reference atom
        #   - CN: known coordination number of complex (default 6)
        # OUTPUT
        #   - nats: list of indices of connected atoms
        ratom = self.getAtom(ind)
        #print('called slow function...')
        # calculates adjacent number of atoms
        nats = []
        for i,atom in enumerate(self.atoms):
            valid = True # flag 
            d = distance(ratom.coords(),atom.coords())
            ## default interatomic radius
            ## for non-metalics
            distance_max = 1.15*(atom.rad+ratom.rad)
            if atom.ismetal() or ratom.ismetal(): 
                if debug:
                    print('metal in  cat ' + str(atom.symbol()) + ' and rat ' +str(ratom.symbol()) )
                ## one the atoms is a metal!
                ## use a longer max for metals
                distance_max = 1.30*(atom.rad+ratom.rad) 
                if d < distance_max and i!=ind:
                    ### trim Hydrogens
                    if atom.symbol() == 'H' or ratom.symbol() == 'H':
                        if debug:
                            print('invalid due to hydrogens: ')
                            print(atom.symbol())
                            print(ratom.symbol())
                        valid = False
                    if d < distance_max and i!=ind and valid:
                        if atom.symbol() == "C":           
                            ## in this case, atom might be intruder C!
                            possible_inds = self.getBondedAtomsnotH(ind) ## bonded to metal
                            if len(possible_inds)>CN:
                                metal_prox = sorted(possible_inds,key=lambda x: self.getDistToMetal(x,ind))
                               
                                allowed_inds = metal_prox[0:CN]
                                if debug:
                                    print('ind: '+str(ind))
                                    print('metal prox:' + str(metal_prox))
                                    print('trimmed to '+str(allowed_inds))
                                    print(allowed_inds)
                                if not i in allowed_inds:
                                    valid = False
                                    if debug:
                                        print('bond rejected based on atom: ' + str(i) + ' not in ' +str(allowed_inds))
                                else:
                                    if debug:
                                        print('Ok based on atom')
                        if ratom.symbol() == "C":  
                            ## in this case, ratom might be intruder C!
                            possible_inds = self.getBondedAtomsnotH(i) ## bonded to metal
                            metal_prox = sorted(possible_inds,key=lambda x: self.getDistToMetal(x,i))
                            if len(possible_inds)>CN:
                                allowed_inds = metal_prox[0:CN]
                                if debug:
                                    print('ind: '+str(ind))
                                    print('metal prox:' + str(metal_prox))
                                    print('trimmed to '+str(allowed_inds))
                                    print(allowed_inds)
                                if not ind in allowed_inds:
                                    valid = False
                                    if debug:
                                        print('bond rejected based on ratom ' + str(ind) + ' with symbol ' + ratom.symbol())
                                else:
                                    if debug:
                                        print('ok based on ratom...')
                else:
                    if debug:
                        print('distance too great')
            if (d < distance_max and i!=ind):
                if valid:
                    if debug:
                        print('Valid atom  ind ' + str(i) + ' (' + atom.symbol() + ') and ' + str(ind) + ' (' + ratom.symbol() + ')')
                        print(' at distance ' + str(d) + ' (which is less than ' + str(distance_max) + ')')
                    nats.append(i)
                else:
                    if debug:
                        print('atom  ind ' + str(i) + ' (' + atom.symbol() + ')')
                        print('has been disallowed from bond with ' + str(ind) + ' (' + ratom.symbol() + ')')
                        print(' at distance ' + str(d) + ' (which would normally be less than ' + str(distance_max) + ')')
                    if d<2 and not atom.symbol() == 'H' and not ratom.symbol() == 'H':
                        print('Error, mol3D could not understand conenctivity in mol' )
        return nats
        
        
    ## Gets atoms bonded to a specific atom using the molecular graph, or creates it
    #
    #  @param self The object pointer
    #  @param ind Index of reference atom
    #  @return List of indices of bonded atoms
    def getBondedAtomsSmart(self,ind):
        if not len(self.graph):
            self.createMolecularGraph()
        return list(np.nonzero(np.ravel(self.graph[ind]))[0])

    ## Gets non-H atoms bonded to a specific atom
    #
    #  Otherwise identical to getBondedAtoms().
    #  @param self The object pointer
    #  @param ind Index of reference atom
    #  @return List of indices of bonded atoms
    def getBondedAtomsnotH(self,ind):
        ratom = self.getAtom(ind)
        # calculates adjacent number of atoms
        nats = []
        for i,atom in enumerate(self.atoms):
            d = distance(ratom.coords(),atom.coords())
            distance_max = 1.15*(atom.rad+ratom.rad)
            if atom.ismetal() or ratom.ismetal(): 
                distance_max = 1.30*(atom.rad+ratom.rad) 
            else:
                distance_max = 1.15*(atom.rad+ratom.rad)
            if (d < distance_max and i!=ind and atom.sym!='H'):
                nats.append(i)
        return nats

    ## Gets atom that is furthest from the molecule COM along a given direction and returns the corresponding distance
    #  @param self The object pointer
    #  @param uP Search direction
    #  @return Distance
    def getfarAtomdir(self,uP):
        dd = 1000.0
        atomc = [0.0,0.0,0.0]
        for atom in self.atoms:
            d0 = distance(atom.coords(),uP)
            if d0 < dd:
                dd = d0
                atomc = atom.coords()
        return distance(self.centermass(),atomc)

    ## Gets H atoms in molecule
    #  @param self The object pointer
    #  @return List of atom3D objects of H atoms
    def getHs(self):
        hlist = []
        for i in range(self.natoms):
            if self.getAtom(i).sym=='H':
                hlist.append(i)
        return hlist

    ## Gets H atoms bonded to specific atom3D in molecule
    #  @param self The object pointer
    #  @param ratom atom3D of reference atom
    #  @return List of atom3D objects of H atoms
    def getHsbyAtom(self,ratom):
        nHs = []
        for i,atom in enumerate(self.atoms):
            if atom.sym == 'H':
                d = distance(ratom.coords(),atom.coords())
                if (d < 1.2*(atom.rad+ratom.rad) and d > 0.01):
                    nHs.append(i)
        return nHs

    ## Gets H atoms bonded to specific atom index in molecule
    #
    #  Trivially equivalent to getHsbyAtom().
    #  @param self The object pointer
    #  @param idx Index of reference atom
    #  @return List of atom3D objects of H atoms
    def getHsbyIndex(self,idx):
        # calculates adjacent number of hydrogens
        nHs = []
        for i,atom in enumerate(self.atoms):
            if atom.sym == 'H':
                d = distance(atom.coords(),self.getAtom(idx).coords())
                if (d < 1.2*(atom.rad+self.getAtom(idx).rad) and d > 0.01):
                    nHs.append(i)
        return nHs

    ## Gets index of closest atom to reference atom.
    #  @param self The object pointer
    #  @param atom0 Index of reference atom
    #  @return Index of closest atom
    def getClosestAtom(self,atom0):
        # INPUT
        #   - atom0: reference atom3D
        # OUTPUT
        #   - idx: index of closest atom to atom0 from molecule
        idx = 0
        cdist = 1000
        for iat,atom in enumerate(self.atoms):
            ds = atom.distance(atom0)
            if (ds < cdist):
                idx = iat
                cdist = ds
        return idx

    ## Gets point that corresponds to mask
    #  @param self The object pointer    
    #  @param mask Identifier for atoms
    #  @return Center of mass of mask
    def getMask(self,mask):
        globs = globalvars()
        elements = globs.elementsbynum()
        # check center of mass
        ats = []
        # loop over entries in mask
        for entry in mask:
            # check for center of mass
            if ('com' in entry.lower()) or ('cm' in entry.lower()):
                return self.centermass()
            # check for range
            elif '-' in entry:
                at0 = entry.split('-')[0]
                at1 = entry.split('-')[-1]
                for i in range(int(at0),int(at1)+1):
                    ats.append(i-1) # python indexing
            elif entry in elements:
                ats += self.findAtomsbySymbol(entry)
            else:
                # try to convert to integer
                try:
                    t = int(entry)
                    ats.append(t-1)
                except:
                    return self.centermass()
        maux = mol3D()
        for at in ats:
            maux.addAtom(self.getAtom(at))
        if maux.natoms==0:
            return self.centermass()
        else:
            return maux.centermass()

    ## Gets index of closest non-H atom to another atom
    #  @param self The object pointer    
    #  @param atom0 atom3D of reference atom
    #  @return Index of closest non-H atom
    def getClosestAtomnoHs(self,atom0):
        idx = 0
        cdist = 1000
        for iat,atom in enumerate(self.atoms):
            ds = atom.distance(atom0)
            if (ds < cdist) and atom.sym!='H':
                idx = iat
                cdist = ds
        return idx
        
    ## Gets distance between two atoms in molecule
    #  @param self The object pointer    
    #  @param idx Index of first atom
    #  @param metalx Index of second atom
    #  @return Distance between atoms
    def getDistToMetal(self,idx,metalx):
        d = self.getAtom(idx).distance(self.getAtom(metalx))
        return d

    ## Gets index of closest non-H atom to another atom
    #
    #  Equivalent to getClosestAtomnoHs() except that the index of the reference atom is specified.
    #  @param self The object pointer    
    #  @param atidx Index of reference atom
    #  @return Index of closest non-H atom
    def getClosestAtomnoHs2(self,atidx):
        idx = 0
        cdist = 1000
        for iat,atom in enumerate(self.atoms):
            ds = atom.distance(self.getAtom(atidx))
            if (ds < cdist) and atom.sym!='H' and iat!=atidx:
                idx = iat
                cdist = ds
        return idx

    ## Initializes OBMol object from a file or SMILES string
    #
    #  Uses the obConversion tool and for files containing 3D coordinates (xyz,mol) and the OBBuilder tool otherwise (smiles).
    #  @param self The object pointer    
    #  @param fst Name of input file
    #  @param convtype Input filetype (xyz,mol,smi)
    #  @param ffclean Flag for FF cleanup of generated structure (default False)
    #  @return OBMol object
    def getOBMol(self,fst,convtype,ffclean=False):
        obConversion = openbabel.OBConversion()
        OBMol = openbabel.OBMol()
        if convtype == 'smistring':
            obConversion.SetInFormat('smi')
            obConversion.ReadString(OBMol,fst)
        else:
            obConversion.SetInFormat(convtype[:-1])
            obConversion.ReadFile(OBMol,fst)
        if 'smi' in convtype:
            OBMol.AddHydrogens()    
            b = openbabel.OBBuilder()
            b.Build(OBMol)
        if ffclean:    
            forcefield = openbabel.OBForceField.FindForceField('mmff94')
            forcefield.Setup(OBMol)
            forcefield.ConjugateGradients(200)
            forcefield.GetCoordinates(OBMol)
        self.OBMol = OBMol
        return OBMol

    ## Removes attributes from mol3D object
    #  @param self The object pointer    
    def initialize(self):
        self.atoms = []
        self.natoms = 0
        self.mass = 0
        self.size = 0
        self.graph = []

    ## Calculates the largest distance between atoms of two molecules
    #  @param self The object pointer    
    #  @param mol mol3D of second molecule
    #  @return Largest distance between atoms in both molecules
    def maxdist(self,mol):
        # INPUT
        #   - mol: second molecule
        # OUTPUT
        #   - maxd: maximum distance between atoms of the 2 molecules
        maxd = 0
        for atom1 in mol.atoms:
            for atom0 in self.atoms:
                if (distance(atom1.coords(),atom0.coords()) > maxd):
                    maxd = distance(atom1.coords(),atom0.coords())
        return maxd

    ## Calculates the smallest distance between atoms of two molecules
    #  @param self The object pointer    
    #  @param mol mol3D of second molecule
    #  @return Smallest distance between atoms in both molecules
    def mindist(self,mol):
        # INPUT
        #   - mol: second molecule
        # OUTPUT
        #   - mind: minimum distance between atoms of the 2 molecules
        mind = 1000
        for atom1 in mol.atoms:
            for atom0 in self.atoms:
                if (distance(atom1.coords(),atom0.coords()) < mind):
                    mind = distance(atom1.coords(),atom0.coords())
        return mind

    ## Calculates the smallest distance between atoms in the molecule
    #  @param self The object pointer    
    #  @return Smallest distance between atoms in molecule
    def mindistmol(self):
        mind = 1000
        for ii,atom1 in enumerate(self.atoms):
            for jj,atom0 in enumerate(self.atoms):
                d = distance(atom1.coords(),atom0.coords())
                if  (d < mind) and ii!=jj:
                    mind = distance(atom1.coords(),atom0.coords())
        return mind

    ## Calculates the smallest distance from atoms in the molecule to a given point
    #  @param self The object pointer   
    #  @param point List of coordinates of reference point 
    #  @return Smallest distance to point
    def mindisttopoint(self,point):
        mind = 1000
        for atom1 in self.atoms:
            d = distance(atom1.coords(),point)
            if (d < mind):
                mind = d
        return mind

    ## Calculates the smallest distance between non-H atoms of two molecules
    # 
    #  Otherwise equivalent to mindist().
    #  @param self The object pointer    
    #  @param mol mol3D of second molecule
    #  @return Smallest distance between non-H atoms in both molecules
    def mindistnonH(self,mol):
        mind = 1000
        for atom1 in mol.atoms:
            for atom0 in self.atoms:
                if (distance(atom1.coords(),atom0.coords()) < mind):
                    if (atom1.sym!='H' and atom0.sym!='H'):
                        mind = distance(atom1.coords(),atom0.coords())
        return mind

    ## Calculates the size of the molecule, as quantified by the max. distance between atoms and the COM.
    #  @param self The object pointer    
    #  @return Molecule size (max. distance between atoms and COM)
    def molsize(self):
        maxd = 0
        cm = self.centermass()
        for atom in self.atoms:
            if distance(cm,atom.coords()) > maxd:
                maxd = distance(cm,atom.coords())
        return maxd

    ## Checks for overlap with another molecule
    # 
    #  Compares pairwise atom distances to 0.85*sum of covalent radii
    #  @param self The object pointer    
    #  @param mol mol3D of second molecule
    #  @param silence Flag for printing warnings
    #  @return Flag for overlap
    def overlapcheck(self,mol,silence):
        overlap = False
        for atom1 in mol.atoms:
            for atom0 in self.atoms:
                if (distance(atom1.coords(),atom0.coords()) < 0.85*(atom1.rad + atom0.rad)):
                    overlap = True
                    if not (silence):
                        print "#############################################################"
                        print "!!!Molecules might be overlapping. Increase distance!!!"
                        print "#############################################################"
                    break
        return overlap

    ## Checks for overlap with another molecule with increased tolerance
    # 
    #  Compares pairwise atom distances to 1*sum of covalent radii
    #  @param self The object pointer    
    #  @param mol mol3D of second molecule
    #  @return Flag for overlap
    def overlapcheckh(self,mol):
        overlap = False
        for atom1 in mol.atoms:
            for atom0 in self.atoms:
                if (distance(atom1.coords(),atom0.coords()) < 1.0):
                    overlap = True
                    break
        return overlap
        
        
    ## Prints xyz coordinates to stdout
    # 
    #  To write to file (more common), use writexyz() instead.
    #  @param self The object pointer
    def printxyz(self):
        for atom in self.atoms:
            xyz = atom.coords()
            ss = "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
            print ss

    ## Load molecule from xyz file
    # 
    #  Consider using getOBMol, which is more general, instead.
    #  @param self The object pointer    
    #  @param filename Filename
    def readfromxyz(self,filename):
        self.graph = [] 
        fname = filename.split('.xyz')[0]
        f = open(fname+'.xyz','r')
        s = f.read().splitlines()
        f.close()
        for line in s[2:]:
            l = filter(None,line.split(None))
            if len(l) > 3:
                atom = atom3D(l[0],[float(l[1]),float(l[2]),float(l[3])])
                self.addAtom(atom)

    ## Computes RMSD between two molecules
    # 
    #  Note that this routine does not perform translations or rotations to align molecules.
    #
    #  To do so, use geometry.kabsch().
    #  @param self The object pointer  
    #  @param mol2 mol3D of second molecule
    #  @return RMSD between molecules, NaN if molecules have different numbers of atoms
    def rmsd(self,mol2):
        Nat0 = self.natoms
        Nat1 = mol2.natoms
        if (Nat0 != Nat1):
            print "ERROR: RMSD can be calculated only for molecules with the same number of atoms.."
            return NaN
        else:
            rmsd = 0
            for atom0,atom1 in zip(self.getAtoms(),mol2.getAtoms()):
                rmsd += (atom0.distance(atom1))**2
            rmsd /= Nat0
            return sqrt(rmsd)

    ## Checks for overlap within the molecule
    # 
    #  Single-molecule version of overlapcheck().
    #  @param self The object pointer      
    #  @param silence Flag for printing warning
    #  @return Flag for overlap
    #  @return Minimum distance between atoms
    def sanitycheck(self,silence):
        overlap = False
        mind = 1000
        for ii,atom1 in enumerate(self.atoms):
            for jj,atom0 in enumerate(self.atoms):
                if ii!=jj and (distance(atom1.coords(),atom0.coords()) < 0.7*(atom1.rad + atom0.rad)):
                    overlap = True
                    if distance(atom1.coords(),atom0.coords()) < mind:
                        mind = distance(atom1.coords(),atom0.coords())
                    if not (silence):
                        print "#############################################################"
                        print "!!!Molecules might be overlapping. Increase distance!!!"
                        print "#############################################################"
                    break
        return overlap,mind

    ## Translate all atoms by given vector.
    #  @param self The object pointer      
    #  @param dxyz Translation vector
    def translate(self,dxyz):
        for atom in self.atoms:
            atom.translate(dxyz)

    ## Writes xyz file in GAMESS format
    #  @param self The object pointer   
    #  @param filename Filename    
    def writegxyz(self,filename):
        ss = '' # initialize returning string
        ss += "Date:"+time.strftime('%m/%d/%Y %H:%M')+", XYZ structure generated by mol3D Class, "+self.globs.PROGRAM+"\nC1\n"
        for atom in self.atoms:
            xyz = atom.coords()
            ss += "%s \t%.1f\t%f\t%f\t%f\n" % (atom.sym,float(atom.atno),xyz[0],xyz[1],xyz[2])
        fname = filename.split('.gxyz')[0]
        f=open(fname+'.gxyz','w')
        f.write(ss)
        f.close()

    ## Writes xyz file
    #
    #  To print to stdout instead, use printxyz().
    #  @param self The object pointer   
    #  @param filename Filename  
    def writexyz(self,filename):
        ss = '' # initialize returning string
        ss += str(self.natoms)+"\n"+time.strftime('%m/%d/%Y %H:%M')+", XYZ structure generated by mol3D Class, "+self.globs.PROGRAM+"\n"
        for atom in self.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        fname = filename.split('.xyz')[0]
        f=open(fname+'.xyz','w')
        f.write(ss)
        f.close()

    ## Writes xyz file for 2 molecules combined
    #
    #  Used when placing binding molecules.
    #  @param self The object pointer   
    #  @param mol mol3D of second molecule    
    #  @param filename Filename  
    def writemxyz(self,mol,filename):
        ss = '' # initialize returning string
        ss += str(self.natoms+mol.natoms)+"\n"+time.strftime('%m/%d/%Y %H:%M')+", XYZ structure generated by mol3D Class, "+self.globs.PROGRAM+"\n"
        for atom in self.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        for atom in mol.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        fname = filename.split('.xyz')[0]
        f=open(fname+'.xyz','w')
        f.write(ss)
        f.close()

    ## Writes xyz file for 2 molecules separated
    #
    #  Used when placing binding molecules for computation of binding energy.
    #  @param self The object pointer   
    #  @param mol mol3D of second molecule    
    #  @param filename Filename  
    def writesepxyz(self,mol,filename):
        ss = '' # initialize returning string
        ss += str(self.natoms)+"\n"+time.strftime('%m/%d/%Y %H:%M')+", XYZ structure generated by mol3D Class, "+self.globs.PROGRAM+"\n"
        for atom in self.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        ss += "--\n"+str(mol.natoms)+"\n\n"
        for atom in mol.atoms:
            xyz = atom.coords()
            ss += "%s \t%f\t%f\t%f\n" % (atom.sym,xyz[0],xyz[1],xyz[2])
        fname = filename.split('.xyz')[0]
        f=open(fname+'.xyz','w')
        f.write(ss)
        f.close()

    ## Print methods
    #  @param self The object pointer   
    #  @return String with methods
    def __repr__(self):
        # OUTPUT
        #   - ss: string with all methods
        # overloaded function
        """ when calls mol3D object without attribute e.g. t """
        ss = "\nClass mol3D has the following methods:\n"
        for method in dir(self):
            if callable(getattr(self, method)):
                ss += method +'\n'
        return ss
