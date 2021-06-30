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
	
	def __init__(self, use_atom_specific_cutoffs=False):
		# Number of atoms part of proteins
		self.natoms = 0 
		# Number of atoms not part of proteins
		self.nhetatms = 0
		# Number of chains
		self.nchains = 0
		# List of atoms part of proteins
		self.atoms = []
		# List of atoms not part of proteins
		self.hetatms = []
		# List of chains
		self.chains = []
		# List of missing atoms
		self.missing_atoms = []
		# List of missing amino acids
		self.missing_aas = []
