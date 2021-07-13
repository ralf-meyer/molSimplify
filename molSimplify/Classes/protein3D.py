#  @file protein3D.py
#  Defines protein3D class and contains useful manipulation/retrieval routines.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

# imports
from math import sqrt
import os
from molSimplify.Classes.AA3D import AA3D
from molSimplify.Classes.atom3D import atom3D
import gzi
from itertools import chain

# no GUI support for now

class protein3D:
    """Holds information about a protein, used to do manipulations.  Reads
    information from structure file (pdb, cif) or is directly built from
    molsimplify.
    
    """
    
    def __init__(self, pdbfile='undef'):
        # Number of amino acids
        self.naas = 0 
        # Number of atoms not part of proteins
        self.nhetatms = 0
        # Number of chains
        self.nchains = 0
        # Dictionary of amino acids
        self.aas = {}
        # Dictionary of all atoms
        self.atoms = {}
        # Dictionary of heteroatoms
        self.hetatms = {}
        # Dictionary of chains
        self.chains = {}
        # Dictionary of missing atoms
        self.missing_atoms = {}
        # List of missing amino acids
        self.missing_aas = []
        # List of possible AA conformations 
        self.conf = []
        # R value
        self.R = -1
        # Rfree value
        self.Rfree = -1
        # PDB file (if applicable)
        self.pdbfile = pdbfile
        # Holder for metals
        self.metals = False
        # Bonds
        self.bonds = {}
    
    def setAAs(self, aas):
        """ Set amino acids of a protein3D class to different amino acids.
        Parameters
        ----------
            aas : dictionary
                Keyed by AA3D amino acids
                Valued by atom3D atoms that are in the amino acids
        """
        self.aas = aas
        self.naas = len(aas.keys())

    def setHetatms(self, hetatms):
        """ Set heteroatoms of a protein3D class to different heteroatoms.
        Parameters
        ----------
            hetatms : dictionary
                Keyed by heteroatoms
                Valued by a list containing the overarching molecule and chain
        """
        self.hetatms = hetatms
        self.nhetatms = len(hetatms.values())

    def setChains(self, chains):
        """ Set chains of a protein3D class to different chains.

        Parameters
        ----------
            chains : dictionary
                Keyed by desired chain IDs.
                Valued by the list of AA3D amino acids in the chain.
        """
        self.chains = chains
        self.nchains = len(chains.keys())

    def setMissingAtoms(self, missing_atoms):
        """ Set missing atoms of a protein3D class to a new dictionary.

        Parameters
        ----------
            missing_atoms : dictionary
                Keyed by amino acid residues of origin
                Valued by missing atoms
        """
        self.missing_atoms = missing_atoms

    def setMissingAAs(self, missing_aas):
        """ Set missing amino acids of a protein3D class to a new list.

        Parameters
        ----------
            missing_aas : list
                List of missing amino acids.
        """
        self.missing_aas = missing_aas
            
    def setConf(self, conf):
        """ Set possible conformations of a protein3D class to a new list.

        Parameters
        ----------
            conf : list
                List of possible conformations for applicable amino acids.
        """
        self.conf = conf
            
    def setR(self, R):
        """ Set R value of protein3D class.

        Parameters
        ----------
                R : float
                        The desired new R value.
        """
        self.R = R
            
    def setRfree(self, Rfree):
        """ Set Rfree value of protein3D class.

        Parameters
        ----------
                Rfree : float
                        The desired new Rfree value.
        """
        self.Rfree = Rfree
            
    def getMissingAtoms(self):
        """ Get missing atoms of a protein3D class.

        """
        return self.missing_atoms.values()
    
    def getMissingAAs(self):
        """ Get missing amino acid residues of a protein3D class.

        """
        return self.missing_aas
    
    def countAAs(self):
        """ Return the number of amino acid residues in a protein3D class.

        """
        return self.naas

    def findAtom(self, sym="X"):
        """
        Find atoms with a specific symbol that are contained in amino acids.
        
        Parameters
        ----------
            sym: str
                element symbol, default as X.

        Returns
        ----------
            inds: list
                a list of atom index with the specified symbol.
        """
        inds = []
        for ii, atom in enumerate(self.aas.values()):
            if atom.symbol() == sym:
                inds.append(ii)
        return inds

    def findHetAtom(self, sym="X"):
        """
        Find heteroatoms with a specific symbol.
        
        Parameters
        ----------
            sym: str
                element symbol, default as X.

        Returns
        ----------
            inds: list
                a list of atom index with the specified symbol.
        """
        inds = []
        for ii, atom in enumerate(self.hetatms.values()):
            if atom.symbol() == sym:
                inds.append(ii)
        return inds

    def findAA(self, three_lc="XAA"):
        """
        Find amino acids with a specific three-letter code.
        
        Parameters
        ----------
            three_lc: str
                three-letter code, default as XAA.

        Returns
        -------
            inds: list
                a list of amino acid indices with the specified symbol.
        """
        inds = []
        for ii, aa in enumerate(self.aas.keys()):
            if aa.three_lc == three_lc:
                inds.append(ii)
        return inds

    def getChain(self, chain_id):
        """ Takes a chain of interest and turns it into its own protein3D class instance.

        Parameters
        ----------
            chain_id : string
                The letter name of the chain of interest

        Returns
        -------
            p : protein3D
                A protein3D instance consisting of just the chain of interest
        """
        p = protein3D()
        p.setChain(chain[chain_id])
        p.setAAs(p.chain.values())
        p.setR(self.R)
        p.setRfree(self.Rfree)
        missing_aas = []
        for aa in self.missing_aas:
            if aa.chain == chain_id:
                missing_aas.append(aa)
        p.setMissingAAs(missing_aas)
        gone_chain = self.missing_atoms
        for aa in gone_chain.keys():
            if aa.chain != chain_id:
                del gone_chain[aa]
        p.setMissingAtoms(gone_chain)
        gone_hets = self.hetatms
        for het in gone_hets.keys():
            if gone_hets[het][1] != chain_id:
                del gone_hets[het]
        p.setHetatms(gone_hets)
        return p

    def getResidue(self, a_id):
        """ Finds the amino acid residue that the atom is contained in.

        Parameters
        ----------
            a_id : int
                the index of the desired atom whose residue we want to find

        Returns
        -------
            aa : AA3D
                the amino acid residue containing the atom
                returns None if there is no amino acid
        """
        for aa in self.aas.keys():
            if (a_id, atoms[a_id]) in self.aas[aa]:
                return aa
        for aa in self.missing_atoms.keys():
            if (a_id, atoms[a_id]) in self.missing_atoms[aa]:
                return aa
        return None # the atom is a heteroatom

    def stripAtoms(self, atoms_stripped):
        """ Removes certain atoms from the protein3D class instance.
        
        Parameters
        ----------
            atoms_stripped : list
                list of atom3D indices that should be removed
        """
        for aa in self.aas.keys():
            for (a_id, atom) in self.aas[aa]:
                if a_id in atoms_stripped:
                    self.aas[aa].remove((a_id, atom))
                    atoms_stripped.remove(a_id)
        for atoms[h_id] in self.hetatms.keys():
            if h_id in atoms_stripped:
                del self.hetatms[atoms[h_id]]   
                atoms_stripped.remove(h_id)

    def stripHetMol(self, hetmol):
        """ Removes all heteroatoms part of the specified heteromolecule from
            the protein3D class instance.

        Parameters
        ----------
            hetmol : str
                String representing the name of a heteromolecule whose
                heteroatoms should be stripped from the protein3D class instance
        """
        for (h_id, hetatm) in self.hetatms.keys():
            if hetmol in self.hetatms[(h_id, hetatm)]:
                del self.hetatms[(h_id, hetatm)] 

    def findMetal(self, transition_metals_only=True):
        """Find metal(s) in a protein3D class.
        Parameters
        ----------
            transition_metals_only : bool, optional
                Only find transition metals. Default is true.
                
        Returns
        -------
            metal_list : list
                List of indices of metal atoms in protein3D.
        """
        if not self.metals:
            metal_list = []
            for i, atom in enumerate(self.hetatms.keys()): # no metals in AAs
                if atom.ismetal(transition_metals_only=transition_metals_only):
                    metal_list.append(i)
            self.metals = metal_list
        return (self.metals)

    def freezeatom(self, atomIdx):
        """Set the freeze attribute to be true for a given atom3D class.

        Parameters
        ----------
            atomIdx : int
                Index for atom to be frozen.
        """

        self.atoms[atomIdx].frozen = True

    def freezeatoms(self, Alist):
        """Set the freeze attribute to be true for a given set of atom3D classes, 
        given their indices. Preserves ordering, starts from largest index.

        Parameters
        ----------
            Alist : list
                List of indices for atom3D instances to remove.
        """

        for h in sorted(Alist, reverse=True):
            self.freezeatom(h)

    def getAtom(self, idx):
        """Get atom with a given index.

        Parameters
        ----------
            idx : int
                Index of desired atom.

        Returns
        -------
            atom : atom3D
                atom3D class for element at given index.

        """

        return self.atoms[idx]

    def getBoundAAs(self, hetatm):
        """Get a list of amino acids bound to a heteroatom, usually a metal.

        Parameters
        ----------
            hetatm : atom3D
                the desired (hetero)atom origin

        Returns
        -------
            bound_aas : list
                list of AA3D instances of amino acids bound to hetatm
        """
        bound_aas = []
        for b in self.bonds[hetatm]:
            if self.getResidue(b) != None:
                bound_aas.append(self.getResidue(b))
        return bound_aas
    
    def readfrompdb(self, text):
        """ Read PDB into a protein3D class instance.

        Parameters
        ----------
            text : str
                String of path to PDB file. Path may be local or global.
                May also be the text of a PDB file from the internet.
        """
        if '.pdb' in  text: # means this is a filename
            self.pdbfile = text
            fname = filename.split('.pdb')[0]
            f = open(fname + '.pdb', r)
            text = f.read()
        # class attributes
        aas = {}
        hetatms = {}
        atoms = {}
        chains = {}
        missing_atoms = {}
        missing_aas = []
        conf = []
        bonds = {}
        f.close()
        # get R and Rfree values
        temp = text.split("R VALUE            (WORKING SET) : ")
        temp = temp[-1]
        temp = temp.split('\nREMARK   3   FREE R VALUE                     : ')
        R = float(temp[0])
        temp = temp[1].split('\n')
        Rfree = float(temp[0])
        # start getting missing amino acids
        text = text.split("M RES C SSSEQI")
        want = text[-1]
        text = text[0].split('\n')
        split = text[-1]
        want = want.split(split)
        for line in want:
            if line == want[-1]:
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            a = AA3D(l[0], l[1], l[2])
            missing_aas.append(a)
        # start getting missing atoms
        text = text.split("M RES CSSEQI  ATOMS")
        want = text[-1]
        text = text[0].split('\n')
        split = text[-1]
        want = want.split(split)
        for line in want:
            if line == want[-1]: 
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            a = AA3D(l[0], l[1], l[2])
            missing_atoms[a] = []
            for atom in l[3:]:
                missing_atoms[a].append(atom3D(Sym=atom))
        # start getting chains - initialize keys of dictionary
        text = text.split('\nSEQRES')
        for line in text:
            if line == text[-1]:
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            if l[2] not in chains.keys():
                chains[l[2]] = [] # this just gets the letter of the chain
        # start getting amino acids
        text = text.split('\nATOM')
        for line in text:
            if line == text[-1]:
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            a = AA3D(l[2], l[3], l[4], float(l[8]))
            if int(l[8]) != 1 and a not in conf:
                conf.append(a)
            if a not in chains[l[2]] and a not in conf:
                chains[l[2]].append(a)
            if a not in aas.keys():
                aas[a] = []
            atom = atom3D(Sym=l[10], xyz=[l[5], l[6], l[7]], Tfactor=l[9],
                          occup=float(l[8]), greek=l[1])
            aas[a].append((int(l[0]), atom)) # terminal Os may be missing
            atoms[int(l[0])] = atom
            bonds.update(a.bonds)
            if a.prev != None:
                bonds[aa.n].add(aa.prev.c)
            if aa.next != None:
                bonds[aa.c].add(aa.next.n)
        # start getting hetatoms
        text = text.split('\nHETATM')
        for line in text:
            if line == text[-1]:
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            hetatm = atom3D(sym=l[-1], xyz = [l[5], l[6], l[7]], Tfactor=l[9],
                            occup=float(l[8]), greek=l[1])
            if hetatm not in hetatms.keys():
                hetatms[hetatm] = [l[2], l[3]] # [cmpd name, chain]
            atoms[int(l[0])] = hetatm
        # deal with conformations
        for i in range(len(conf)-1):
            if conf[i].chain == conf[i+1].chain and conf[i].id == conf[i+1].id:
                if conf[i].occup >= conf[i+1].occup:
                    chains[conf[i].chain].append(conf[i])
                    # pick chain with higher occupancy or the A chain if tie
                else:
                    chains[conf[i+1].chain].append(conf[i+1])
        # get extra connections
        text = text.split('\nCONECT')
        for line in text:
            if line == text[-1]:
                text = line
                line = line.split('\n')
                line = line[:1]
                text.replace(line, '')
            l = line.split()
            if atoms[int(l[0])] not in bonds.keys():
                bonds[atoms[int(l[0])]] = set()
            for i in l[1:]:
                bonds[atoms[int(l[0])]].add(atoms[int(i)])
        self.setChains(chains)
        self.setAAs(aas)
        self.setHetatms(hetatms)
        self.setMissingAtoms(missing_atoms)
        self.setMissingAAs(missing_aas)
        self.setConf(conf)
        self.setR(R)
        self.setRfree(Rfree)

    def fetch_pdb(self, pdbCode):
        """ API query to fetch a pdb and write it as a protein3D class instance

        Parameters
        ----------
            pdbCode : not sure what this variable is oops
        """
        try:
            import urllib.request as urllib
        except ImportError:
            import urllib                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        

        remoteCode = pdbCode.upper()
        if not os.path.exists(pdb_dir):
            os.mkdir(pdb_dir)
        try:
            filename = urllib.urlretrieve(
                'http://www.rcsb.org/pdb/cgi/export.cgi/' + remoteCode +
                '.pdb.gz?format=PDB&pdbId=' + remoteCode +
                '&compression=gz')[0]
            temp_p = protein3D(pdbfile=filename)
        except:
            print("warning: %s not found.\n"%pdbCode)
        else:
            if (os.path.getsize(filename) > 0): # If 0, pdb code was invalid
                try:
                    abort = 0
                    temp_p.readfrompdb(gzip.open(filename).read())
                    print("fetched: %s"%(pdbCode))
                except IOError:
                    print('aborted')
            else:
                print("warning: %s not valid.\n"%pdbCode)

    def setBonds(self, bonds):
        """Sets the bonded atoms in the protein.
        This is effectively the molecular graph.
         
        Parameters
        ----------
            bonds : dictionary
                Keyed by atom3D atoms in the protein
                Valued by a set consisting of bonded atoms
        """
        self.bonds = bonds

