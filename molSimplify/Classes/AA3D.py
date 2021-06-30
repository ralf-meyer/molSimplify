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
	"""Holds information about an amino acid, ussed to do manipulations.  Reads information from structure file (pdb, cif) or is directly built from molsimplify.
	
	"""
	
	def __init__(self, use_atom_specific_cutoffs=False):
		# List of atom3D objects
		self.atoms = []
		# Number of atoms
		self.natoms = 0
		# 3-letter code of amino acid (in all caps)
		self.name = ""
