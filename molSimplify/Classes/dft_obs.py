x# Written by JP Janet for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## Defines class of DFT observations################
########     will be used to postprocess #################
########   DFT results by measuring ligand properties ####
##########################################################


from molSimplify.Classes.ligand import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.atom3D import *
from molSimplify.Informatics.autocorrelation import *
class dft_observation:
	def __init__(self,name,geopath):
		self.name = name
		self.descriptor = dict()
		self.mol = False
		self.health = False
		self.commnets = list()
		self.geopath = geopath
	def sety(self,y_value):
		self.yvalue = y_value
	def obtain_mol3d(self):
		new_mol = mol3D()
		this_mol = mol3D()
		this_mol.readfromxyz(self.geopath)
		if this_mol.natoms > 0:
			self.mol = this_mol
			self.natoms = this_mol.natoms
			self.health = True
		else:
			self.comments.append('geo file appears empty')
			self.health = False
	def get_descriptor_vector(self,loud):
		
		results_dictionary = generate_all_ligand_misc(self.mol,loud)
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_ax'],'misc','ax')
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_eq'],'misc','eq')
		results_dictionary = generate_all_ligand_autocorrelations(self.mol,depth=3,loud=loud,name=name)
		
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_ax_full'],'f','ax')
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_eq_full'],'f','eq')
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_ax_con'],'c','ax')
		self.append_descriptors(results_dictionary['colnames'],results_dictionary['result_eq_con'],'c','eq')
	def append_descriptors(self,list_of_names,list_of_props,prefix,suffix):
		for names in list_of_names:
			if hasattr(names, '__iter__'):
				names = ["-".join([prefix,str(i),suffix]) for i in names]
				self.descriptor_names += names
			else:
				names = "-".join([prefix,str(names),suffix])
				self.descriptor_names.append(names)
			
			
 		for values in list_of_props:
			if hasattr(values, '__iter__'):
				self.descriptors.extend(values)
			else:
				self.descriptors.append(values)
