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
		
	def identify(self):
		""" States whether the amino acid is (positively/negatively) charged, polar, amphipathic, or hydrophobic.
		
		Returns
		-------
		aa_type : string
			Positively charged, Negatively charged, Polar, Hydrophobic
			
		"""
		if self.name == "ARG" or self.name == "LYS":  return "Positively charged"
		elif self.name == "ASP" or self.name == "GLU":  return "Negatively charged"
		elif self.name in {"GLN", "ASN", "HIS", "SER", "THR", "TYR", "CYS"}: return "Polar"
		else: return "Hydrophobic"
