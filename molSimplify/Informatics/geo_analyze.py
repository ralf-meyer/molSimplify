# Written JP Janet
# for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## Defines class for postprocessing    #############
########     geometries for octahedral    ################
################      TM complexs          ###############
##########################################################
from molSimplify.Classes.ligand import *


def getOctBondDistances(mol):
    ## This function gets
    ## ax and equitorial 
    ## min and max bond lengths
    liglist,ligdents,ligcons = ligand_breakdown(mol)
    ax_ligand_list,eq_ligand_list,ax_natoms_list,eq_natoms_list,ax_con_int_list,eq_con_int_list,ax_con_list,eq_con_list,built_ligand_list=ligand_assign(mol,
                                                                                                                                                        liglist,
                                                                                                                                                        ligdents,
                                                                                                                                                        ligcons,
                                                                                                                                                        False,'g')
    ax_dist = list()
    eq_dist = list()
    for ax_ligs in ax_con_list:
        tempList = list()
        for conatms in ax_ligs:
            tempList.append(distance(mol.getAtom(mol.findMetal()[0]).coords(),mol.getAtom(conatms).coords()))
        ax_dist.append(tempList)
    for eq_ligs in eq_con_list:
        tempList = list()
        for conatms in eq_ligs:
            tempList.append(distance(mol.getAtom(mol.findMetal()[0]).coords(),mol.getAtom(conatms).coords()))
        eq_dist.append(tempList)
    print(ax_dist)
    print(eq_dist)
    return ax_dist, eq_dist
def maximum_ML_dist(mol):

    core = mol.getAtom(mol.findMetal()[0]).coords()
    max_dist = 0
    for atom_inds in mol.getBondedAtomsSmart(mol.findMetal()[0]):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist > max_dist):
            max_dist = dist
    return max_dist
def minimum_ML_dist(mol):
    
    core = mol.getAtom(mol.findMetal()[0]).coords()
    min_dist = 1000
    for atom_inds in mol.getBondedAtomsSmart(mol.findMetal()[0]):
        dist = distance(core,mol.getAtom(atom_inds).coords())
        if (dist < min_dist) and (dist > 0):
            min_dist = dist
    return min_dist
def mean_ML_dist(mol):
    core = mol.getAtom(mol.findMetal()[0]).coords()
    mean_dist = 0.0
    for atom_inds in mol.getBondedAtomsSmart(mol.findMetal()[0]):
	dist = distance(core,mol.getAtom(atom_inds).coords())
        mean_dist += float(dist)
    mean_dist = mean_dist/float(len(mol.getBondedAtomsSmart(mol.findMetal()[0])))
    return mean_dist


