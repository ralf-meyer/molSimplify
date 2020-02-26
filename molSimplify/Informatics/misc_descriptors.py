import glob
import string
import sys
import os
import numpy as np
import math
import random
import string
import numpy
from molSimplify.Scripts.geometry import *
#from molSimplify.Scripts.nn_prep import *
from molSimplify.Classes.atom3D import *
from molSimplify.Classes.globalvars import globalvars
from molSimplify.Classes.mol3D import*
from molSimplify.Classes.ligand import *
from molSimplify.Classes.ligand import ligand_breakdown#, ligand_assign
from molSimplify.Classes.ligand import ligand_assign_consistent as ligand_assign
from molSimplify.Informatics.graph_analyze import *


def generate_all_ligand_misc(mol, loud, custom_ligand_dict=False, force_legacy=False):
    # custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
    ##                                    ax_con_int_list ,eq_con_int_list
    # with types: eq/ax_ligand_list list of mol3D
    # eq/ax_con_int_list list of list/tuple of int e.g,  [[1,2] [1,2]]
    # use force_legacy to get kier indices/MCDL
    result_ax = list()
    result_eq = list()
    if force_legacy:
        colnames = ['dent', 'maxDEN', 'ki', 'tki', 'charge']
    else:
        colnames = ['dent', 'charge']
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
            mol, liglist, ligdents, ligcons, loud, name=False)
    else:
        ax_ligand_list = custom_ligand_dict["ax_ligand_list"]
        eq_ligand_list = custom_ligand_dict["eq_ligand_list"]
        ax_con_int_list = custom_ligand_dict["ax_con_int_list"]
        eq_con_int_list = custom_ligand_dict["eq_con_int_list"]
    # count ligands
    n_ax = len(ax_ligand_list)
    n_eq = len(eq_ligand_list)
    # allocate
    result_ax_dent = False
    result_eq_dent = False
    result_ax_maxdelen = False
    result_eq_maxdelen = False
    result_ax_ki = False
    result_eq_ki = False
    result_ax_tki = False
    result_eq_tki = False
    result_ax_charge = False
    result_eq_charge = False
    # loop over axial ligands
    if n_ax > 0:
        for i in range(0, n_ax):
            ax_ligand_list[i].mol.convert2OBMol()
            if not (i == 0):
                result_ax_dent += ax_ligand_list[i].dent
                if force_legacy:
                    result_ax_maxdelen += get_lig_EN(
                        ax_ligand_list[i].mol, ax_con_int_list[i])
                    result_ax_ki += kier(ax_ligand_list[i].mol)
                    result_ax_tki += get_truncated_kier(
                        ax_ligand_list[i].mol, ax_con_int_list[i])
                result_ax_charge += ax_ligand_list[i].mol.OBMol.GetTotalCharge()
            else:
                result_ax_dent = ax_ligand_list[i].dent
                if force_legacy:
                    result_ax_maxdelen = get_lig_EN(
                        ax_ligand_list[i].mol, ax_con_int_list[i])
                    result_ax_ki = kier(ax_ligand_list[i].mol)
                    result_ax_tki = get_truncated_kier(
                        ax_ligand_list[i].mol, ax_con_int_list[i])
                result_ax_charge = ax_ligand_list[i].mol.OBMol.GetTotalCharge()
        # average axial results
        result_ax_dent = np.divide(result_ax_dent, n_ax)
        if force_legacy:
            result_ax_maxdelen = np.divide(result_ax_maxdelen, n_ax)
            result_ax_ki = np.divide(result_ax_ki, n_ax)
            result_ax_tki = np.divide(result_ax_tki, n_ax)
        result_ax_charge = np.divide(result_ax_charge, n_ax)

    # loop over eq ligands
    if n_eq > 0:
        for i in range(0, n_eq):
            eq_ligand_list[i].mol.convert2OBMol()
            if not (i == 0):
                result_eq_dent += eq_ligand_list[i].dent
                if force_legacy:
                    result_eq_maxdelen += get_lig_EN(
                        eq_ligand_list[i].mol, eq_con_int_list[i])
                    result_eq_ki += kier(eq_ligand_list[i].mol)
                    result_eq_tki += get_truncated_kier(
                        eq_ligand_list[i].mol, eq_con_int_list[i])
                result_eq_charge += eq_ligand_list[i].mol.OBMol.GetTotalCharge()
            else:
                result_eq_dent = eq_ligand_list[i].dent
                if force_legacy:
                    result_eq_maxdelen = get_lig_EN(
                        eq_ligand_list[i].mol, eq_con_int_list[i])
                    result_eq_ki = kier(eq_ligand_list[i].mol)
                    result_eq_tki = get_truncated_kier(
                        eq_ligand_list[i].mol, eq_con_int_list[i])
                result_eq_charge = eq_ligand_list[i].mol.OBMol.GetTotalCharge()
        # average eq results
        result_eq_dent = np.divide(result_eq_dent, n_eq)
        if force_legacy:
            result_eq_maxdelen = np.divide(result_eq_maxdelen, n_eq)
            result_eq_ki = np.divide(result_eq_ki, n_eq)
            result_eq_tki = np.divide(result_eq_tki, n_eq)
        result_eq_charge = np.divide(result_eq_charge, n_eq)
        # save the results
    result_ax.append(result_ax_dent)
    if force_legacy:
        result_ax.append(result_ax_maxdelen)
        result_ax.append(result_ax_ki)
        result_ax.append(result_ax_tki)
    result_ax.append(result_ax_charge)

    result_eq.append(result_eq_dent)
    if force_legacy:
        result_eq.append(result_eq_maxdelen)
        result_eq.append(result_eq_ki)
        result_eq.append(result_eq_tki)
    result_eq.append(result_eq_charge)

    results_dictionary = {'colnames': colnames,
                          'result_ax': result_ax, 'result_eq': result_eq}
    return results_dictionary

def get_lig_EN(mol, connection_atoms):
    # calculate the maximum abs electronegativity
    # difference between connection atom an all
    # neighbors
    max_EN = 0
    globs = globalvars()
    for atoms in connection_atoms:
        this_atoms_neighbors = mol.getBondedAtomsSmart(atoms)
        for bound_atoms in this_atoms_neighbors:
            this_EN = float(globs.endict()[mol.getAtom(atoms).symbol(
            )]) - float(globs.endict()[mol.getAtom(bound_atoms).symbol()])
            if (abs(this_EN) >= max_EN):
                max_EN = this_EN
    return max_EN
