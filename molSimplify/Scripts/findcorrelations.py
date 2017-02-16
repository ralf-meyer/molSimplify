# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT
import os, sys, copy
import glob, re, math, random, string, numpy, pybel
from math import pi
from molSimplify.Scripts.geometry import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.mol3D import*
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Informatics.graph_analyze import *

def get_descriptors():
	# returns allowed descriptors
	avail_descriptors  = ['kier','truncated_kier','randic','truncated_randic',
				'con_atom_type','con_atom_en','max_con_atom_delta_en',
				'autocorrelation']
	return(avail_descriptors)
def check_descriptors(args):
	clean_descriptors = list()
	allowed_descriptors = get_descriptors()
	if args.descriptors: ## check if the descriptors are 
			     ## correctly given
		if args.debug:	
			print('descriptors are' + str(args.descriptors))
		for descriptors in args.descriptors:
			if descriptors in allowed_descriptors:
				clean_descriptors.append(descriptors)
			else:
				print('ignoring invalid descriptor ' + str(descriptors))
				print('allowed values are  '+ str(allowed_descriptors))
	else:
		clean_descriptors = allowed_descriptors
	n_desc = len(clean_descriptors)
	if args.debug:
		print('using ' + str(n_desc) + ' descriptors: ' + str(clean_descriptors) )
	if n_desc > 0:
		return clean_descriptors
	else:
		print('Error, no valid descriptors! using default instead')
		return(allowed_descriptors)
def analysis_supervisor(args):
	## main routine to control
	## analysis process
	descriptors = check_descriptors(args)
	print('more to come!')


