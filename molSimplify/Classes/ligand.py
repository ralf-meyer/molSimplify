# @file ligand.py
#  Defines ligand class for postprocessing DFT results by measuring ligand properties
#
#  Written by JP Janet for HJK Group
#
#  Dpt of Chemical Engineering, MIT

from molSimplify.Classes.mol3D import *


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

    # map betweem
    # int and ext indcies
    # Obtain the ligand from the complex mol3D object
    # @param self The object pointer
    def obtain_mol3d(self):
        this_mol = mol3D()
        this_ext_int_dict = dict()
        j = 0
        # the old routine where all atoms in the master_mol are gone through from 0 to natoms-1
        # for i in range(0, self.master_mol.natoms):
        #     if i in self.index_list:
        # the new rountine where the indices are taken out directly. This way the order of atoms is preserved
        for i in self.index_list:
            this_mol.addAtom(self.master_mol.getAtom(i))
            this_ext_int_dict.update({i: j})
            j += 1  # keep count of how many are added
        self.mol = this_mol
        self.ext_int_dict = this_ext_int_dict

    # Truncate ligand about connecting atoms
    #  @param self The object pointer
    #  @param con_atoms The connection atom indices
    #  @param hops Number of bonds to truncate after
    #  @return Truncated mol3D object
    def obtain_truncation(self, con_atoms, hops):
        self.trunc_mol = mol3D()
        added_list = list()
        for connections in con_atoms:
            hopped = 0
            active_set = [connections]
            while hopped < hops:
                hopped += 1
                new_active_set = list()
                for this_atom in active_set:
                    this_atoms_neighbors = self.master_mol.getBondedAtomsSmart(
                        this_atom)
                    for bound_atoms in this_atoms_neighbors:
                        if (bound_atoms in self.index_list) and (bound_atoms not in added_list):
                            self.trunc_mol.addAtom(
                                self.master_mol.getAtomSmart(bound_atoms))
                            added_list.append(bound_atoms)
                    [new_active_set.append(element)
                     for element in this_atoms_neighbors]
                active_set = new_active_set
