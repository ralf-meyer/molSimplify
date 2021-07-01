#  @file protein3D.py
#  Defines protein3D class and contains useful manipulation/retrieval routines.
#
#  Written by HJK Group
#
#  Dpt of Chemical Engineering, MIT

# imports
from math import sqrt
import os

# no GUI support for now

class protein3D:
	"""Holds information about a protein, used to do manipulations.  Reads information from structure file (pdb, cif) or is directly built from molsimplify.
	
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
		# List of atoms not part of proteins
		self.hetatms = {}
		# List of chains
		self.chains = []
		# Dictionary of missing atoms
		self.missing_atoms = {}
		# List of missing amino acids
		self.missing_aas = []
		self.pdbfile = pdbfile
	
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
				Keyed by the molecule that the heteroatoms together form
				Valued by heteroatoms
		"""
		self.hetatms = hetatms
		self.nhetatms = len(hetatms.values())

	def setChains(self, chains):
		""" Set chains of a protein3D class to different chains.
		Parameters
		----------
			chains : list
				List of desired chain IDs.
		"""
		self.chains = chains
		self.nchains = len(chains)

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

	def readfrompdb(self, filename):
		""" Read PDB into a protein3D class instance.
		Parameters
		----------
			filename : string
				String of path to PDB file. Path may be local or global.
		
		Returns
		-------
			p : protein3D
				A protein3D instance created from the provided PDB file.
		"""
		self.pdbfile = filename
		fname = filename.split('.pdb')[0]
		f = open(fname + '.pdb', r)
		text = f.read()
		p = protein3D()
		# class attributes
		aas = {}
		hetatms = {}
		chains = set() # because we only want distinct chains
		missing_atoms = {}
		missing_aas = []
		f.close()
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
		text = text.split("M RES C SSSEQI  ATOMS")
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
		# start getting chains
		text = text.split('\nSEQRES')
		for line in text:
			if line == text[-1]:
				text = line
				line = line.split('\n')
				line = line[:1]
				text.replace(line, '')
			l = line.split()
			chains.add(l[2]) # this just gets the letter of the chain
			chains = list(chains)
		# start getting amino acids
		text = text.split('\nATOM')
		for line in text:
			if line == text[-1]:
				text = line
				line = line.split('\n')
				line = line[:1]
				text.replace(line, '')
			l = line.split()
			a = AA3D(l[2], l[4], l[5])
			if a not in aas.keys(): aas[a] = []
			atom = atom3D(Sym=l[-1], xyz=[l[5], l[6], l[7]])
			aas[a].append(atom)
		# start getting hetatoms
		text = text.split('\nHETATM')
		for line in text:
			if line == text[-1]:
				text = line
				line = line.split('\n')
				line = line[:1]
				text.replace(line, '')
			l = line.split()
			if l[2] not in hetatms.keys(): hetatms[l[2]] = [] # l[2] is the name of a compound
			hetatm = atom3D(sym=l[-1], xyz = [l[5], l[6], l[7]])
			hetatms[l[2]].append(hetatm)
		p.setChains(chains)
		p.setAAs(aas)
		p.setHetatms(hetatms)
		p.setMissingAtoms(missing_atoms)
		p.setMissingAAs(missing_aas)
		return p


