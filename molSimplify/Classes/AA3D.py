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
		
		
	def identify(self):
		""" States whether the amino acid is (positively/negatively) charged, polar, or hydrophobic.
		
		Returns
		-------
		aa_type : string
			Positively charged, Negatively charged, Polar, Hydrophobic
			
		"""
		if self.name == "ARG" or self.name == "LYS":  return "Positively charged"
		elif self.name == "ASP" or self.name == "GLU":  return "Negatively charged"
		elif self.name in {"GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS"}: return "Polar"
		else: return "Hydrophobic"
		
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
		
	
