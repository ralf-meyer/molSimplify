# @file ligand.py
#  Defines ligand class for postprocessing DFT results by measuring ligand properties
#
#  Written by JP Janet for HJK Group
#
#  Dpt of Chemical Engineering, MIT



# Ligand class for postprocessing DFT results by measuring ligand properties
class ligand:
    # Constructor
    #  @param self The object pointer
    #  @param master_mol A mol3D complex to extract ligand from
    #  @param index_list A list of indices of the ligand atoms
    #  @paran dent The denticity of the ligand
    def __init__(self, master_mol, index_list, dent):
        self.master_mol = master_mol
        self.index_list = index_list
        self.dent = dent
        self.ext_int_dict = dict()  # store


    ### Deprecated.
    # map betweem
    # int and ext indcies
    # Obtain the ligand from the complex mol3D object
    # @param self The object pointer
    # def obtain_mol3d(self):
    #     from molSimplify.Classes.mol3D import mol3D
    #     this_mol = mol3D()
    #     this_ext_int_dict = dict()
    #     j = 0
    #     # the old routine where all atoms in the master_mol are gone through from 0 to natoms-1
    #     # for i in range(0, self.master_mol.natoms):
    #     #     if i in self.index_list:
    #     # the new rountine where the indices are taken out directly. This way the order of atoms is preserved
    #     for i in self.index_list:
    #         this_mol.addAtom(self.master_mol.getAtom(i))
    #         this_ext_int_dict.update({i: j})
    #         j += 1  # keep count of how many are added
    #     self.mol = this_mol
    #     self.ext_int_dict = this_ext_int_dict
