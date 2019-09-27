# @file ligand.py
#  Defines ligand class for postprocessing DFT results by measuring ligand properties
#
#  Written by JP Janet for HJK Group
#
#  Dpt of Chemical Engineering, MIT
import numpy as np
from itertools import combinations
from molSimplify.Classes.mol3D import mol3D


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

    ## Deprecated.
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


##########

# Extract axial and equitorial components of a octahedral complex
#  @param mol The mol3D object for the complex
#  @param liglist List of ligands
#  @return ligdents List of ligand dents
#  @return ligcons List of ligand connection indices (in mol)
def ligand_breakdown(mol, flag_loose=False, BondedOct=False, silent=True):
    # this function takes an octahedral
    # complex and returns ligands
    loud = False
    metal_index = mol.findMetal()[0]
    bondedatoms = mol.getBondedAtomsSmart(metal_index, oct=True)
    # print('!!!!!boundatoms', bondedatoms)
    #	print('from get oct' + str(bondedatoms))
    #	print('***\n')
    bonded_atom_symbols = [mol.getAtom(i).symbol() for i in bondedatoms]
    if not silent:
        print('result of ligand ligand_breakdown', bonded_atom_symbols)
    counter = 0
    liglist = []
    ligdents = []
    ligcons = []
    for atom in bondedatoms:
        # print('this atom type is ' + mol.getAtom(atom).symbol())
        # print('conection number ' + str(atom) + " of " + str(bondedatoms))
        fragment = mol.findsubMol(atom, metal_index)
        this_cons = [x for x in fragment if (x in bondedatoms)]
        # print('fragment',fragment)
        # print('this_cons',this_cons)
        unique = True
        for i, unique_ligands in enumerate(liglist):
            if sorted(fragment) == sorted(unique_ligands):
                unique = False
                matched = i
        if unique:
            liglist.append(fragment)
            ligdents.append(1)
            ligcons.append(this_cons)
        else:
            ligdents[matched] += 1
    return liglist, ligdents, ligcons


def ligand_assign(mol, liglist, ligdents, ligcons, loud=False, name=False):
    valid = True
    # loud = False
    pentadentate = False
    metal_index = mol.findMetal()[0]
    built_ligand_list = list()
    lig_natoms_list = list()
    unique_ligands = list()
    ligand_counts = list()
    all_ligand_counts = [0, 0, 0, 0, 0, 0]
    ligand_records = list()
    ax_con_int_list = list()
    eq_con_int_list = list()
    ax_natoms_list = list()
    eq_natoms_list = list()
    n_ligs = len(liglist)
    max_dent = max(ligdents)
    min_dent = min(ligdents)
    if loud:
        print('********************************************')
        print(("n_ligs = " + str(n_ligs)))
        print(("max d = " + str(max_dent)))
        print(("min_dent = " + str(min_dent)))
        print(("ligand list is" + str(liglist)))
        print(('denticities are  ' + str(ligdents)))
    if (max(ligdents) == 4) and (min(ligdents) != 1):
        valid = False
        print(('bad denticities: ' + str(ligdents)))
        print(('min denticities: ' + str(min(ligdents))))
    if max(ligdents) > 4:
        #### Handling of pentadentate ligands goes here. #####
        if max(ligdents) == 5 and min(ligdents) == 1:
            pentadentate = True
        elif max(ligdents) == 6 and min(ligdents) == 6:
            hexadentate = True
        else:
            valid = False
            print('bad denticities: ' + str(ligdents))
            print('max denticities: ' + str(min(ligdents)))
    if n_ligs > 3 and min(ligdents) > 1:
        valid = False
        print(('too many ligs ' + str((n_ligs))))
    eq_lig_list = list()
    ax_lig_list = list()
    ax_con_list = list()
    eq_con_list = list()
    for i, ligand_indices in enumerate(liglist):
        this_ligand = ligand(mol, ligand_indices, ligdents[i])
        this_ligand.obtain_mol3d()
        # lig_natoms_list.append(this_ligand.mol.natoms) ## old one with obtain_mol3d
        built_ligand_list.append(this_ligand)
        lig_natoms_list.append(len(this_ligand.index_list))
    for j, built_ligs in enumerate(built_ligand_list):
        # test if ligand is unique
        sl = [atom.symbol() for atom in built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
        # _sl = [atom.symbol() for atom in built_ligs.mol.getAtoms()] ## old one with obtain_mol3d
        if loud:
            print(('checking lig ' + str(j) + ' : ' + str(sl)))
        unique = 1
        for i, other_sl in enumerate(unique_ligands):
            if sorted(sl) == sorted(other_sl):
                # duplicate
                unique = 0
                ligand_counts[i] += 1
        if unique == 1:
            unique_ligands.append(sl)
            ligand_counts.append(1)
            ligand_records.append(j)
    # loop to bin ligands:
    for j, built_ligs in enumerate(built_ligand_list):
        # test if ligand is unique
        sl = [atom.symbol() for atom in built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
        unique = 1
        for i, other_sl in enumerate(unique_ligands):
            if sorted(sl) == sorted(other_sl):
                # duplicate
                # print(i,ligand_counts[i])
                all_ligand_counts[j] = ligand_counts[i]

    if loud:
        print(('unique ligands' + str(unique_ligands)))
        print(('ligand counts' + str(ligand_counts)))
        print(('ligand records ' + str(ligand_records)))
        print((str(max(ligand_counts)) +
               ' is the max and min in  ' + str(min(ligand_counts))))
    n_unique_ligs = len(unique_ligands)
    if (n_ligs == 3) or (n_ligs == 4):  # most common case,
        # one/two equitorial and 2 axial mono
        # or three bidentate
        for i, ligs in enumerate(liglist):
            if ligdents[i] == 1 and min_dent == 1:  # anything with equitorial monos will
                # have higher than 4 n_ligs
                ax_lig_list.append(i)
                if loud:
                    print(('choosing ' + str(i) + ' as ax based on dent =1'))
                ax_con_list.append(ligcons[i])
            if (ligdents[i] >= 2) and (min_dent == 1):
                eq_lig_list.append(i)
                if loud:
                    print(('choosing lig ' + str(i) + ' as eq based on high dent'))
                eq_con_list.append(ligcons[i])
        if (n_ligs == 3) and (min_dent == max_dent):
            if n_unique_ligs == 1:
                # take any 2, they are all the same
                if loud:
                    print('triple bidentate case')
                ax_lig_list.append(0)
                eq_lig_list.append(1)
                eq_lig_list.append(2)
                ax_con_list.append(ligcons[0])
                eq_con_list.append(ligcons[1])
                eq_con_list.append(ligcons[2])
            elif min_dent == 2 and max_dent == 2 and n_ligs == 3 and not n_unique_ligs == 1:
                # this is a hetero/bidentate case
                for i, ligs in enumerate(liglist):
                    if all_ligand_counts[i] == 2:
                        eq_lig_list.append(i)
                        eq_con_list.append(ligcons[i])
                    elif all_ligand_counts[i] == 1:
                        ax_lig_list.append(i)
                        ax_con_list.append(ligcons[i])
    elif (n_ligs == 6):  # all mono  case,
        minz = 500
        maxz = -500
        if loud:
            print('monodentate case')
        allowed = list(range(0, 6))
        not_eq = list()
        for j, built_ligs in enumerate(built_ligand_list):
            this_z = sum([mol.getAtom(ii).coords()[2]
                          for ii in ligcons[j]]) / len(ligcons[j])
            if this_z < minz:
                minz = this_z
                bot_lig = j
                bot_con = ligcons[j]
            if loud:
                print(('updating bot axial to ' + str(bot_lig)))
            if this_z > maxz:
                maxz = this_z
                top_lig = j
                top_con = ligcons[j]
            if loud:
                print(('updating top axial to ' + str(top_lig)))
        not_eq.append(bot_lig)
        not_eq.append(top_lig)

        allowed = [x for x in allowed if ((x not in not_eq))]
        if len(allowed) != 4:
            print(('error in decomp of monodentate case!', allowed))
        eq_lig_list = allowed
        eq_con_list = [ligcons[i] for i in allowed]
        ax_lig_list = [top_lig, bot_lig]
        ax_con_list = [top_con, bot_con]
        if loud:
            print(('geometric eq_list ' + str(eq_lig_list)))
            print(('geometric ax_list ' + str(eq_lig_list)))
        if (max(ligand_counts) != 4) or (min(ligand_counts) != 2):
            if loud:
                print('not a 4-6 case')
            if (max(ligand_counts) == 6):
                if loud:
                    print('6-homoleptic, using geo values')
            # ax=ligand_records[ligand_counts.index(6)]
            # eq_lig=ligand_records[ligand_counts.index(6)]
            else:
                if loud:
                    print('monodentates not the same, using geo values ')
                    print(ligand_counts)
                    print(unique_ligands)
        elif n_unique_ligs == 2:
            if loud:
                print('this is a  4-6 case')
            allowed = list(range(0, 6))
            ax_lig_list = [i for i in allowed if (all_ligand_counts[i] == 2)]
            eq_lig_list = [i for i in allowed if (all_ligand_counts[i] == 4)]
            ax_con_list = [ligcons[i] for i in ax_lig_list]
            eq_con_list = [ligcons[i] for i in eq_lig_list]
    elif n_ligs == 2 and pentadentate:
        #### Handling for pentadentate scaffolds ####
        minz = 500
        maxz = -500
        if loud:
            print('pentadentate case')
        allowed = [0, 1]
        not_eq = list()
        for j, built_ligs in enumerate(built_ligand_list):
            if len(ligcons[j]) == 1:
                #### This is the axial ligand ####
                # print(j, 'axial lig')
                top_lig = j
                top_con = ligcons[j]
                not_eq.append(top_lig)
            else:
                pentadentate_coord_list = np.array([mol.getAtom(
                    ii).coords() for ii in ligcons[j]])
                ##### Adjusting this so that by default, any 4 within the same plane will be assigned as eq. ###
                if loud:
                    print('pentadentate coord LIST!')
                    print(pentadentate_coord_list)
                point_combos = combinations([0, 1, 2, 3, 4], 4)
                error_list = []
                combo_list = []
                for i, combo in enumerate(point_combos):
                    combo_list.append(list(combo))
                    A = []
                    b = []
                    for point_num in combo:
                        coordlist = pentadentate_coord_list[point_num]
                        A.append([coordlist[0], coordlist[1], 1])
                        b.append(coordlist[2])
                    ##### This code builds the best fit plane between 4 points,
                    ##### Then calculates the variance of the 4 points with respect to the plane
                    ##### The 4 that have the least variance are flagged as the eq plane.
                    mat_b = np.matrix(b).T
                    mat_A = np.matrix(A)
                    try:
                        fit = (mat_A.T * mat_A).I * mat_A.T * mat_b
                        errors = np.squeeze(np.array(mat_b - mat_A * fit))
                        error_var = np.var(errors)
                        error_list.append(error_var)
                    except:
                        error_list.append(0)
                if loud:
                    print('combos below')
                    print(combo_list)
                    print('errors next, argmin combo selected')
                    print(error_list)
                not_ax_points = combo_list[np.argmin(error_list)]
                if len(set(not_ax_points)) != 4:
                    print('The equatorial plane is not being assigned correctly. Please check.')
                    sardines
                else:
                    bot_idx = list(set(range(5)) - set(not_ax_points))[0]
                    if loud:
                        print('This is bot_idx', bot_idx)
                    bot_lig = j
                    bot_con = [ligcons[j][bot_idx]]

        allowed = list(set(allowed) - set(not_eq))
        if loud:
            print('this is the allowed list', allowed, not_eq)
        eq_lig_list = allowed
        eq_con_list = [
            list(set([ligcons[i] for i in allowed][0]) - set(top_con) - set(bot_con))]
        ax_lig_list = [top_lig, bot_lig]
        ax_con_list = [top_con, bot_con]
        if loud:
            print('con lists', eq_con_list, ax_con_list)
        ###########################################################################################
        # In the above, the pentadentate ligand is classified as both axial and equatorial.       #
        # The lc atoms are decided by the z-position. Thus the pentadentate ligand has 4 eq-lc    #
        # and 1 ax-lc. Currently should be able to check this and set that up.                    #
        ###########################################################################################
    elif n_ligs == 1 and hexadentate:
        allowed = [0, 1]
        if loud:
            print('hexadentate case')
        not_eq = list()
        for j, built_ligs in enumerate(built_ligand_list):
            hexadentate_coord_list = np.array([mol.getAtom(
                ii).coords() for ii in ligcons[j]])
            ##### Adjusting this so that by default, any 4 within the same plane will be assigned as eq. ###
            if loud:
                print('hexadentate coord LIST!')
                print(hexadentate_coord_list)
            # point_combos = combinations([1,2,3,4,5,6],4)
            pair_combos = list(combinations([0, 1, 2, 3, 4, 5], 2))
            angle_list = []
            pair_list = []
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(hexadentate_coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(hexadentate_coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            argsort_angle_list = np.squeeze(np.array(angle_list)).argsort()[-3:][::-1]
            point_combos = [pair_list[argsort_angle_list[0]] + pair_list[argsort_angle_list[1]],
                            pair_list[argsort_angle_list[1]] + pair_list[argsort_angle_list[2]],
                            pair_list[argsort_angle_list[2]] + pair_list[argsort_angle_list[0]]]
            error_list = []
            combo_list = []
            fitlist = []
            for i, combo in enumerate(point_combos):
                combo_list.append(combo)
                A = []
                b = []
                for point_num in combo:
                    coordlist = hexadentate_coord_list[point_num]
                    A.append([coordlist[0], coordlist[1], 1])
                    b.append(coordlist[2])
                ##### This code builds the best fit plane between 4 points,
                ##### Then calculates the variance of the 4 points with respect to the plane
                ##### The 4 that have the least variance are flagged as the eq plane.
                mat_b = np.matrix(b).T
                mat_A = np.matrix(A)
                fit = (mat_A.T * mat_A).I * mat_A.T * mat_b
                fitlist.append(fit)
                errors = np.squeeze(np.array(mat_b - mat_A * fit))
                error_var = np.var(errors)
                error_list.append(error_var)
            if loud:
                print('combos below')
                print(combo_list)
                print('errors next, argmin combo selected')
                print(error_list)
            best_fit_planes = np.squeeze(np.array(error_list)).argsort()[:3]
            perpdist = []
            perpcombo = []
            for fitnum, best_fit in enumerate(best_fit_planes):
                temp_fit = fitlist[best_fit]
                temp_combo = combo_list[best_fit]
                perpcombo.append(int(best_fit))
                temp_ax = set(range(0, 6)) - set(temp_combo)
                ax_dist = []
                for point_num in temp_ax:
                    coordlist = hexadentate_coord_list[point_num]
                    planez = [coordlist[0], coordlist[1], 1] * temp_fit
                    # planez = temp_fit[0] * coordlist[0] + temp_fit[1] * coordlist[1] + fit[2]
                    plane_coords = [coordlist[0], coordlist[1], planez]
                    adjusted_coords = [coordlist[0], coordlist[1], coordlist[2]]
                    squared_dist = np.sum((np.array(adjusted_coords) - np.array(plane_coords)) ** 2)
                    dist = np.squeeze(np.sqrt(squared_dist))
                    if loud:
                        print('dist', dist)
                    ax_dist.append(dist)
                perpdist.append(np.mean(ax_dist))
            if loud:
                print("Perpendicular distance is", perpdist, perpcombo, len(perpdist), len(best_fit_planes))
            not_ax_points = combo_list[perpcombo[np.argmax(np.array(perpdist))]]
            if len(set(not_ax_points)) != 4:
                print('The equatorial plane is not being assigned correctly. Please check.')
                sardines
            else:
                bot_idx = list(set(range(6)) - set(not_ax_points))[0]
                top_idx = list(set(range(6)) - set(not_ax_points))[1]
                if loud:
                    print('This is bot_idx', bot_idx)
                bot_lig = j
                top_lig = j
                bot_con = [ligcons[j][bot_idx]]
                top_con = [ligcons[j][top_idx]]
        allowed = list(set(allowed) - set(not_eq))
        if loud:
            print('this is the allowed list', allowed, not_eq)
        eq_lig_list = [top_lig]
        eq_con_list = eq_con_list = [
            list(set([ligcons[0][i] for i in not_ax_points]))]
        ax_lig_list = [top_lig, bot_lig]
        ax_con_list = [top_con, bot_con]
        if loud:
            print('con lists', eq_con_list, ax_con_list)

    ############### DONE WITH CLASSIFICATION ######
    # ax_lig=ligand_records[ligand_counts.index(2)]
    # eq_lig=ligand_records[ligand_counts.index(4)]
    ax_ligand_list = [built_ligand_list[i] for i in ax_lig_list]
    eq_ligand_list = [built_ligand_list[i] for i in eq_lig_list]
    if loud and valid:
        print(('lig_nat_list', lig_natoms_list))
        print(('eq_liq is ind ', eq_lig_list))
        print(('ax_liq is ind ', ax_lig_list))
        print(('ax built lig [0] ext ind :' +
               str(list(built_ligand_list[ax_lig_list[0]].ext_int_dict.keys()))))
        if len(ax_lig_list) > 1:
            print(('ax built lig [1] ext ind :' +
                   str(list(built_ligand_list[ax_lig_list[1]].ext_int_dict.keys()))))
        print(('eq built lig [0] ext ind: ' +
               str(list(built_ligand_list[eq_lig_list[0]].ext_int_dict.keys()))))
        print(('eq_con is ' + str((eq_con_list))))
        print(('ax_con is ' + str((ax_con_list))))
    if name:  ## TODO: this part might be broken. Maybe not of much use though.
        for i, ax_ligand in enumerate(ax_ligand_list):
            if not os.path.isdir('ligands'):
                os.mkdir('ligands')
            ax_ligand.mol.writexyz('ligands/' + name +
                                   '_' + str(i) + '_ax.xyz')
        for i, eq_ligand in enumerate(eq_ligand_list):
            if not os.path.isdir('ligands'):
                os.mkdir('ligands')
            eq_ligand.mol.writexyz('ligands/' + name +
                                   '_' + str(i) + '_eq.xyz')
    for j, ax_con in enumerate(ax_con_list):
        current_ligand_index_list = built_ligand_list[ax_lig_list[j]].index_list
        ax_con_int_list.append([current_ligand_index_list.index(i) for i in ax_con])
        # ax_con_int_list.append(
        #     [built_ligand_list[ax_lig_list[j]].ext_int_dict[i] for i in ax_con])  # convert to interal index ## old one with obtain_mol3d
    for j, eq_con in enumerate(eq_con_list):
        current_ligand_index_list = built_ligand_list[eq_lig_list[j]].index_list
        eq_con_int_list.append([current_ligand_index_list.index(i) for i in eq_con])
        # eq_con_int_list.append(
        #     [built_ligand_list[eq_lig_list[j]].ext_int_dict[i] for i in eq_con])  # convert to interal index ## old one with obtain_mol3d
    if loud:
        print(('int eq ' + str(eq_con_int_list)))
        print(('ext eq ' + str(eq_con_list)))
        print('**********************************************')
    for ax_lig in ax_lig_list:
        ax_natoms_list.append(lig_natoms_list[ax_lig])
    for eq_lig in eq_lig_list:
        eq_natoms_list.append(lig_natoms_list[eq_lig])
    return ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, built_ligand_list


def ligand_assign_consistent(mol, liglist, ligdents, ligcons, loud=False, name=False, use_z = False):
    ####### This ligand assignment code handles octahedral complexes consistently.
    ####### It should be able to assign any octahedral complex
    angle_cutoff = 130
    valid = True
    hexadentate = False
    pentadentate = False
    metal_index = mol.findMetal()[0]
    built_ligand_list = list()
    lig_natoms_list = list()
    unique_ligands = list()
    ligand_counts = list()
    all_ligand_counts = [0, 0, 0, 0, 0, 0]
    ligand_records = list()
    ax_con_int_list = list()
    eq_con_int_list = list()
    ax_natoms_list = list()
    eq_natoms_list = list()
    n_ligs = len(liglist)
    max_dent = max(ligdents)
    min_dent = min(ligdents)
    if loud:
        print('********************************************')
        print(("n_ligs = " + str(n_ligs)))
        print(("max d = " + str(max_dent)))
        print(("min_dent = " + str(min_dent)))
        print(("ligand list is" + str(liglist)))
        print(('denticities are  ' + str(ligdents)))
    if max(ligdents) > 4:
        #### Handling of pentadentate ligands goes here. #####
        if max(ligdents) == 5 and min(ligdents) == 1:
            pentadentate = True
        elif max(ligdents) == 6 and min(ligdents) == 6:
            hexadentate = True
        else:
            valid = False
            print('bad denticities: ' + str(ligdents))
            print('max denticities: ' + str(min(ligdents)))
    if n_ligs > 3 and min(ligdents) > 1:
        valid = False
        print(('too many ligs ' + str((n_ligs))))
    eq_lig_list = list()
    ax_lig_list = list()
    ax_con_list = list()
    eq_con_list = list()
    symbol_list = list()
    for i, ligand_indices in enumerate(liglist):
        this_ligand = ligand(mol, ligand_indices, ligdents[i])
        this_ligand.obtain_mol3d()
        built_ligand_list.append(this_ligand)
        lig_natoms_list.append(len(this_ligand.index_list))
    for j, built_ligs in enumerate(built_ligand_list):
        # test if ligand is unique
        sl = [atom.symbol() for atom in built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
        # _sl = [atom.symbol() for atom in built_ligs.mol.getAtoms()] ## old one with obtain_mol3d
        symbol_list.append(sorted(sl))
        if loud:
            print(('checking lig ' + str(j) + ' : ' + str(sl)))
        unique = 1
        for i, other_sl in enumerate(unique_ligands):
            if sorted(sl) == sorted(other_sl):
                # duplicate
                unique = 0
                ligand_counts[i] += 1
        if unique == 1:
            unique_ligands.append(sl)
            ligand_counts.append(1)
            ligand_records.append(j)
    # loop to bin ligands:
    for j, built_ligs in enumerate(built_ligand_list):
        # test if ligand is unique
        sl = [atom.symbol() for atom in built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
        unique = 1
        for i, other_sl in enumerate(unique_ligands):
            print('othersl', other_sl)
            if sorted(sl) == sorted(other_sl):
                # duplicate
                # print(i,ligand_counts[i])
                all_ligand_counts[j] = ligand_counts[i]
    ### Here, obtain a flat last of all connecting atoms (should be of length 6 for octahedral complexes)
    flat_ligcons = [item for sublist in ligcons for item in sublist]
    ### Obtain coordinates for the connecting atoms. Flat coord list ends up being used for comparisons.
    flat_coord_list = np.array([mol.getAtom(ii).coords() for ii in flat_ligcons])
    if loud:
        print('consistent coord LIST!')
        print(flat_coord_list)
    ### Below, take all combinations of two atoms, and measure their angles through the metal
    pair_combos = list(combinations([0, 1, 2, 3, 4, 5], 2))
    angle_list = []
    pair_list = []
    for i, pair in enumerate(pair_combos):
        pair_list.append(list(pair))
        p1 = np.squeeze(np.array(flat_coord_list[list(pair)[0]]))
        p2 = np.squeeze(np.array(flat_coord_list[list(pair)[1]]))
        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
        if loud:
            print('pair of atoms, then angle', pair, angle)
        angle_list.append(angle)
    argsort_angle_list = np.squeeze(np.array(angle_list)).argsort()[-3:][::-1]
    ### Get the 3 largest angles through the metal, which should define the x, y, and z planes of symmetry
    ### Then define those planes of symmetry with the atom indices (point_combos)
    point_combos = [pair_list[argsort_angle_list[0]] + pair_list[argsort_angle_list[1]],
                    pair_list[argsort_angle_list[1]] + pair_list[argsort_angle_list[2]],
                    pair_list[argsort_angle_list[2]] + pair_list[argsort_angle_list[0]]]
    #### Next, fit each of these planes with a best fit plane. 
    #### The plane that fits best will be the equatorial plane by default. 
    #### In some cases, the plane must be overruled (i.e seesaw tetradentates)
    #### For special cases like the seesaw, there is consistent handling such that
    #### there is consistent behavior for the ones within the same ligand class.
    error_list = []
    combo_list = []
    fitlist = []
    for i, combo in enumerate(point_combos):
        combo_list.append(combo)
        if loud:
            print('combo',combo)
        A = []
        b = []
        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
        temp_BL = []
        for point_num in combo:
            coordlist = flat_coord_list[point_num]
            A.append([coordlist[0], coordlist[1], 1])
            b.append(coordlist[2])
        mat_b = np.matrix(b).T
        mat_A = np.matrix(A)
        try:
            fit = (mat_A.T * mat_A).I * mat_A.T * mat_b
            fitlist.append(fit)
            errors = np.squeeze(np.array(mat_b - mat_A * fit))
            error_var = np.var(errors)
            error_list.append(error_var)
        except:
            error_list.append(0) ### perfect fit plane may suffer matrix singularity issues.
    if loud:
        print('combos below')
        print(combo_list)
        print('errors next, argmin combo selected')
        print(error_list)
    #### Find the points that correspond to the best fit plane through 4 points.
    #### Eq points are used later. Eq points has the number of the connection atoms 
    #### across from each other. It pulls from a list of 0, 1, 2, 3, 4, 5, and gets 
    #### which 4 are the atoms that are in the equatorial plane.
    eq_points = combo_list[np.argmin(np.array(error_list))]
    if loud:
        print(('unique ligands' + str(unique_ligands)))
        print(('ligand counts' + str(ligand_counts)))
        print(('ligand records ' + str(ligand_records)))
        print((str(max(ligand_counts)) +
               ' is the max and min in  ' + str(min(ligand_counts))))
    n_unique_ligs = len(unique_ligands)
    if (n_ligs == 6):  # all mono  case,
        allowed = range(0, 6)
        print('monodentate case')
        if n_unique_ligs == 1:
            print('homoleptic monodentate')
            eq_lig_list = range(0, 4)
            eq_con_list = [ligcons[j] for j in eq_lig_list]
            ax_lig_list = range(4, 6)
            ax_con_list = [ligcons[j] for j in ax_lig_list]
        elif n_unique_ligs == 2:
            print('monodentate 5+1')
            if loud:
                print(ligand_counts,unique_ligands)
            eq_lig_list = list()
            if use_z:
                minz = 500
                maxz = -500
                if loud:
                    print('monodentate case')
                allowed = list(range(0, 6))
                not_eq = list()
                for j, built_ligs in enumerate(built_ligand_list):
                    this_z = sum([mol.getAtom(ii).coords()[2]
                                  for ii in ligcons[j]]) / len(ligcons[j])
                    if this_z < minz:
                        minz = this_z
                        bot_lig = j
                        bot_con = ligcons[j]
                    if loud:
                        print(('updating bot axial to ' + str(bot_lig)))
                    if this_z > maxz:
                        maxz = this_z
                        top_lig = j
                        top_con = ligcons[j]
                    if loud:
                        print(('updating top axial to ' + str(top_lig)))
                not_eq.append(bot_lig)
                not_eq.append(top_lig)
                allowed = [x for x in allowed if ((x not in not_eq))]
                if len(allowed) != 4:
                    print(('error in decomp of monodentate case!', allowed))
                eq_lig_list = allowed
                eq_con_list = [ligcons[i] for i in allowed]
                ax_lig_list = [top_lig, bot_lig]
                ax_con_list = [top_con, bot_con]
            else:
                if max(ligand_counts) == 5:
                    for i, ligand_count in enumerate(ligand_counts):
                        temp_unique = unique_ligands[i]
                        for j, built_ligs in enumerate(built_ligand_list):
                            sym_list = [atom.symbol() for atom in
                                        built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
                            if sym_list != temp_unique:
                                continue
                            elif (ligand_count == 5) and len(eq_lig_list) < 4:
                                #### In the 5+1 monodentate case, 4 of the 5 are assigned to be equatorial, 
                                #### last two are assigned to be axial
                                eq_lig_list.append(j)
                                eq_con_list.append(ligcons[j])
                            elif len(ax_lig_list) < 2:
                                ax_lig_list.append(j)
                                ax_con_list.append(ligcons[j])
                else: ### This can be either seesaw configuration or planar configuration, 4+2 cases
                    four_repeats = list()
                    for i, ligand_count in enumerate(ligand_counts):
                        temp_unique = unique_ligands[i]
                        for j, built_ligs in enumerate(built_ligand_list):
                            sym_list = [atom.symbol() for atom in
                                        built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
                            if sym_list != temp_unique:
                                continue
                            elif (ligand_count == 4) and len(four_repeats) < 4:
                                four_repeats.append(j)
                    if loud:
                        print('this is four repeats',four_repeats)
                    four_repeats_cons = [ligcons[j] for j in four_repeats]
                    pair_combos = list(combinations([0, 1, 2, 3], 2))
                    angle_list = []
                    pair_list = []
                    coord_list = np.array([mol.getAtom(ii[0]).coords() for ii in four_repeats_cons])
                    for i, pair in enumerate(pair_combos):
                        pair_list.append(list(pair))
                        p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                        p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        # print('v1v2',v1u,v2u)
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        if loud:
                            print('pair of atoms, then angle', pair, angle)
                        angle_list.append(angle)
                    #### Seesaws will have only 1 ~180 degree angle, whereas planar ligands will have two. 
                    #### Thus, after measuring the angles for all tetradentate connecting atoms through the metal,
                    #### looking at the angle of the first (not zeroeth) element tells us whether or not we have a seesaw.
                    test_angle = np.sort(np.array(angle_list))[::-1][1]
                    if test_angle < angle_cutoff:
                        seesaw = True
                        #### In the seesaw, the two furthest apart are denoted as axial. 
                        #### The equatorial plane consists of 2 seesaw connecting atoms, and two monodentates.
                        axial_pair = [four_repeats_cons[val] for val in list(pair_list[np.argmax(np.array(angle_list))])]
                    else:
                        seesaw = False
                        #### In the planar set, the two monodentates are axial, and tetradentate is equatorial.
                        axial_pair = list(set(allowed) - set(four_repeats))
                    if not seesaw:
                        print('all repeated monodentates in eq plane')
                        eq_lig_list = four_repeats
                        eq_con_list = four_repeats_cons  ### ADDED
                        ax_lig_list = axial_pair
                        ax_con_list = [ligcons[axial_pair[0]], ligcons[axial_pair[1]]]
                    else:  # 2 points in eq plane, seesaw case
                        ### currently this is inconsistent with 
                        print('two repeated monodentates in eq plane')
                        eq_plane_cons = list(set(flat_ligcons) - set([val[0] for val in axial_pair]))
                        eq_con_list = [ligcons[j] for j in range(len(ligcons)) if ligcons[j] in eq_plane_cons]
                        eq_lig_list = [i for i in range(len(ligcons)) if ligcons[i] in eq_plane_cons]
                        ax_lig_list = [i for i in range(len(ligcons)) if ligcons[i] not in eq_plane_cons]
                        ax_con_list = [ligcons[j] for j in range(len(ligcons)) if ligcons[j] not in eq_plane_cons]
        elif n_unique_ligs == 3:
            print('monodentate 4+1+1')
            #### Need to identify if in seesaw-style configuration or planar configuration
            if use_z:
                minz = 500
                maxz = -500
                if loud:
                    print('monodentate case')
                allowed = list(range(0, 6))
                not_eq = list()
                for j, built_ligs in enumerate(built_ligand_list):
                    this_z = sum([mol.getAtom(ii).coords()[2]
                                  for ii in ligcons[j]]) / len(ligcons[j])
                    if this_z < minz:
                        minz = this_z
                        bot_lig = j
                        bot_con = ligcons[j]
                    if loud:
                        print(('updating bot axial to ' + str(bot_lig)))
                    if this_z > maxz:
                        maxz = this_z
                        top_lig = j
                        top_con = ligcons[j]
                    if loud:
                        print(('updating top axial to ' + str(top_lig)))
                not_eq.append(bot_lig)
                not_eq.append(top_lig)
                allowed = [x for x in allowed if ((x not in not_eq))]
                if len(allowed) != 4:
                    print(('error in decomp of monodentate case!', allowed))
                eq_lig_list = allowed
                eq_con_list = [ligcons[i] for i in allowed]
                ax_lig_list = [top_lig, bot_lig]
                ax_con_list = [top_con, bot_con]
            else: 
                if max(ligand_counts) == 4:
                    four_repeats = list()
                    for i, ligand_count in enumerate(ligand_counts):
                        temp_unique = unique_ligands[i]
                        for j, built_ligs in enumerate(built_ligand_list):
                            sym_list = [atom.symbol() for atom in
                                        built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
                            if sym_list != temp_unique:
                                continue
                            elif (ligand_count == 4) and len(four_repeats) < 4:
                                four_repeats.append(j)
                    if loud:
                        print('this is four repeats',four_repeats)
                    four_repeats_cons = [ligcons[j] for j in four_repeats]
                    pair_combos = list(combinations([0, 1, 2, 3], 2))
                    angle_list = []
                    pair_list = []
                    coord_list = np.array([mol.getAtom(ii[0]).coords() for ii in four_repeats_cons])
                    for i, pair in enumerate(pair_combos):
                        pair_list.append(list(pair))
                        p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                        p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        # print('v1v2',v1u,v2u)
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        if loud:
                            print('pair of atoms, then angle', pair, angle)
                        angle_list.append(angle)
                    #### Seesaws will have only 1 ~180 degree angle, whereas planar ligands will have two. 
                    #### Thus, after measuring the angles for all tetradentate connecting atoms through the metal,
                    #### looking at the angle of the first (not zeroeth) element tells us whether or not we have a seesaw.
                    test_angle = np.sort(np.array(angle_list))[::-1][1]
                    if test_angle < angle_cutoff:
                        seesaw = True
                        #### In the seesaw, the two furthest apart are denoted as axial. 
                        #### The equatorial plane consists of 2 seesaw connecting atoms, and two monodentates.
                        axial_pair = [four_repeats_cons[val] for val in list(pair_list[np.argmax(np.array(angle_list))])]
                    else:
                        seesaw = False
                        #### In the planar set, the two monodentates are axial, and tetradentate is equatorial.
                        axial_pair = list(set(allowed) - set(four_repeats))
                    if not seesaw:
                        print('all repeated monodentates in eq plane')
                        eq_lig_list = four_repeats
                        eq_con_list = four_repeats_cons  ### ADDED
                        ax_lig_list = axial_pair
                        ax_con_list = [ligcons[axial_pair[0]], ligcons[axial_pair[1]]]
                    else:  # 2 points in eq plane, seesaw case
                        ### currently this is inconsistent with 
                        print('two repeated monodentates in eq plane')
                        eq_plane_cons = list(set(flat_ligcons) - set([val[0] for val in axial_pair]))
                        eq_con_list = [ligcons[j] for j in range(len(ligcons)) if ligcons[j] in eq_plane_cons]
                        eq_lig_list = [i for i in range(len(ligcons)) if ligcons[i] in eq_plane_cons]
                        ax_lig_list = [i for i in range(len(ligcons)) if ligcons[i] not in eq_plane_cons]
                        ax_con_list = [ligcons[j] for j in range(len(ligcons)) if ligcons[j] not in eq_plane_cons]
                else:
                    three_repeats = list()
                    two_repeats = list()
                    for i, ligand_count in enumerate(ligand_counts):
                        temp_unique = unique_ligands[i]
                        for j, built_ligs in enumerate(built_ligand_list):
                            sym_list = [atom.symbol() for atom in
                                        built_ligs.master_mol.getAtomwithinds(built_ligs.index_list)]
                            if sym_list != temp_unique:
                                continue
                            elif (ligand_count == 3) and len(three_repeats) < 3:
                                three_repeats.append(j)
                            elif (ligand_count == 2) and len(two_repeats) < 2:
                                two_repeats.append(j)
                    if loud:
                        print('this is three repeats',three_repeats)
                    three_repeats_cons = [ligcons[j] for j in three_repeats]
                    pair_combos = list(combinations([0, 1, 2], 2))
                    angle_list = []
                    pair_list = []
                    coord_list = np.array([mol.getAtom(ii[0]).coords() for ii in three_repeats_cons])
                    for i, pair in enumerate(pair_combos):
                        pair_list.append(list(pair))
                        p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                        p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        # print('v1v2',v1u,v2u)
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        if loud:
                            print('pair of atoms, then angle', pair, angle)
                        angle_list.append(angle)
                    test_angle = np.sort(np.array(angle_list))[::-1][0]
                    if test_angle < angle_cutoff:
                        planar = False
                        if len(two_repeats) == 2: ### will only be true for 3+2+1 case
                            print('now in the 3+2+1 case finding the opposite...')
                            two_repeats_cons = [ligcons[j] for j in two_repeats]
                            coord_list_two = np.array([mol.getAtom(ii[0]).coords() for ii in two_repeats_cons])
                            coord_list_three = np.array([mol.getAtom(ii[0]).coords() for ii in three_repeats_cons])
                            angle_list = []
                            tri_bi_pair = []
                            m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                            pair_combos = [(0, 0),(0, 1),(0, 2),(1, 0),(1, 1),(1, 2)]
                            for k, pair in enumerate(pair_combos):
                                p1 = np.squeeze(np.array(coord_list_two[pair[0]]))
                                p2 = np.squeeze(np.array(coord_list_three[pair[1]]))
                                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                                angle_list.append(angle)

                            #### If not planar, one connection atom (for the monodentate ligand) will be opposite
                            #### the connections for the tridentate ligands. Assigning this atom as the opposite monodentate.
                            #### If not planar, considering the plane with the bidentate to be equatorial.
                            top_2_angles = np.squeeze(np.array(angle_list)).argsort()[-2:][::-1]
                            first_pair = pair_combos[top_2_angles[0]]
                            second_pair = pair_combos[top_2_angles[1]]
                            conlist = list()
                            # conlist.append(three_repeats_cons[+two_repeats_cons[tri_bi_pair[top_2_angles[0]][1]]
                            conlist += two_repeats_cons[first_pair[0]]+three_repeats_cons[first_pair[1]]
                            conlist += two_repeats_cons[second_pair[0]]+three_repeats_cons[second_pair[1]]
                            if loud:
                                print('eq points reassigned',conlist)
                            eq_points = [j for j in allowed if ligcons[j][0] in conlist]
                    else:
                        planar = True
                        mono_dentate_idx_set = list(set(allowed)-set(three_repeats))
                        coord_list = np.array([mol.getAtom(ii[0]).coords() for ii in three_repeats_cons])
                        m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                        mono_con_list = []
                        for mono_dentate_idx in mono_dentate_idx_set:
                            current_con = ligcons[mono_dentate_idx]
                            mono_con_list.append(current_con[0])
                            p1 = np.array(mol.getAtom(current_con[0]).coords())
                            angle_list = []
                            for k, tri_con in enumerate(three_repeats_cons):
                                p2 = np.squeeze(np.array(coord_list[k]))
                                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                                angle_list.append(angle)
                            test_angle = np.sort(np.array(angle_list))[::-1][0]
                            #### If the tridentate is planar, only one of the bidentate connecting atoms
                            #### will be across the tridentate. Thus the equatorial plane will be defined
                            #### as the planar tridentate plus one connection atom for the bidentate
                            #### The monodentate and the second bidentate connecting atoms are denoted axial.
                            if test_angle > angle_cutoff:
                                monodentate_eq_con = current_con
                                monodentate_eq_idx = mono_dentate_idx
                    if planar:
                        eq_lig_list = three_repeats+[monodentate_eq_idx]
                        eq_con_list = three_repeats_cons+[monodentate_eq_con]  ### ADDED
                        ax_con_list = [[val] for val in mono_con_list if val not in monodentate_eq_con]
                        # ax_lig_list = list(set(mono_dentate_idx_set)-set([monodentate_eq_idx]))
                        ax_lig_list = [i for i in range(len(ligcons)) if ligcons[i] not in eq_con_list]
                    else:
                        ### any equatorial plane will have 2 of the 3 equivalent ones, so take the eq plane
                        eq_lig_list = eq_points
                        eq_con_list = [ligcons[j] for j in eq_lig_list]
                        # ax_lig_list = list(set(allowed)-set(eq_points))
                        ax_lig_list = [val for i, val in enumerate(allowed) if val not in eq_points]
                        ax_con_list = [ligcons[j] for j in ax_lig_list]
        else:  ### with more than 4 ligands, the ax/eq breaks down. Just going to stick with def of eq plane.
            print('monodentate more than 3 unique')
            eq_ligcons = list(set([flat_ligcons[j] for j in eq_points]))
            eq_lig_list = eq_points
            eq_con_list = [ligcons[j] for j in eq_lig_list]
            ax_lig_list = list(set(allowed)-set(eq_points))
            ax_con_list = [ligcons[j] for j in ax_lig_list]
    elif (n_ligs == 4): # 2+2+1+1 and 3+1+1+1 cases
        allowed = range(0,4)
        if max(ligdents) == 3:
            print('3+1+1+1 case')
            tridentate_ligand_idx = np.argmax(ligdents)
            tridentate_cons = ligcons[tridentate_ligand_idx]
            mono_dentate_idx_set = list(set(range(len(ligdents)))-set([tridentate_ligand_idx]))
            monodentate_cons = [ligcons[val] for val in mono_dentate_idx_set]
            pair_combos = list(combinations([0, 1, 2], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            #### The tridentate can either be planar or not. If planar, it will have at least one
            #### angle that is close to 180 between connecting atoms, through the metal.
            test_angle = np.sort(np.array(angle_list))[::-1][0]
            if test_angle < angle_cutoff:
                planar = False
                full_angle_list = []
                full_con_list = []
                tri_full_con_list = []
                for mono_dentate_idx in mono_dentate_idx_set:
                    current_con = ligcons[mono_dentate_idx]
                    p1 = np.array(mol.getAtom(current_con[0]).coords())
                    if loud:
                        print('monocon',current_con,tridentate_cons)
                    angle_list = []
                    for k, tri_con in enumerate(tridentate_cons):
                        p2 = np.squeeze(np.array(coord_list[k]))
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        angle_list.append(angle)
                    test_angle = np.sort(np.array(angle_list))[::-1][0]
                    #### If the tridentate is planar, only one of the bidentate connecting atoms
                    #### will be across the tridentate. Thus the equatorial plane will be defined
                    #### as the planar tridentate plus one connection atom for the bidentate
                    #### The monodentate and the second bidentate connecting atoms are denoted axial.
                    if test_angle > angle_cutoff:
                        full_angle_list.append(test_angle)
                        full_con_list.append(current_con[0])
                        tri_full_con_list.append(tridentate_cons[np.squeeze(np.array(full_angle_list)).argsort()[-3:][::-1][0]])
                opposite_angle_indices = np.squeeze(np.array(full_angle_list)).argsort()[-3:][::-1]
                mono_lig1 = mono_dentate_idx_set[opposite_angle_indices[0]]
                mono_lig2 = mono_dentate_idx_set[opposite_angle_indices[1]]
                tri_lig_con = [tri_full_con_list[opposite_angle_indices[0]],tri_full_con_list[opposite_angle_indices[1]]]
                tri_lig_ax_con = [val for val in ligcons[tridentate_ligand_idx] if val not in tri_lig_con]
            else:
                print('planar')
                planar = True
                coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons])
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                mono_con_list = []
                for mono_dentate_idx in mono_dentate_idx_set:
                    print(mono_dentate_idx,mono_dentate_idx_set)
                    current_con = ligcons[mono_dentate_idx]
                    mono_con_list.append(current_con[0])
                    p1 = np.array(mol.getAtom(current_con[0]).coords())
                    angle_list = []
                    for k, tri_con in enumerate(tridentate_cons):
                        p2 = np.squeeze(np.array(coord_list[k]))
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        angle_list.append(angle)
                    test_angle = np.sort(np.array(angle_list))[::-1][0]
                    #### If the tridentate is planar, only one of the bidentate connecting atoms
                    #### will be across the tridentate. Thus the equatorial plane will be defined
                    #### as the planar tridentate plus one connection atom for the bidentate
                    #### The monodentate and the second bidentate connecting atoms are denoted axial.
                    if test_angle > angle_cutoff:
                        monodentate_eq_con = current_con
                        monodentate_eq_idx = mono_dentate_idx
            if planar:
                eq_lig_list = [tridentate_ligand_idx, monodentate_eq_idx]
                eq_con_list = [tridentate_cons, monodentate_eq_con]  ### ADDED
                ax_lig_list = [val for val in mono_dentate_idx_set if val not in [monodentate_eq_idx]]
                ax_con_list = [[val] for val in mono_con_list if val not in monodentate_eq_con]
            else:
                ### any equatorial plane will have 2 of the 3 equivalent ones, so take the eq plane
                ax_mono_lig = list(set(allowed)-set([tridentate_ligand_idx,mono_lig1,mono_lig2]))[0]
                eq_lig_list = [tridentate_ligand_idx, mono_lig1, mono_lig2]
                eq_con_list = [tri_lig_con,ligcons[mono_lig1],ligcons[mono_lig2]]
                ax_lig_list = [tridentate_ligand_idx, ax_mono_lig]
                ax_con_list = [tri_lig_ax_con, ligcons[ax_mono_lig]]
        elif max(ligdents) == 2:
            #### Need to handle case with both equatorial and triple-bidentate style.
            print('2+2+1+1 case')
            allowed = range(0, 4)
            bidentate_ligand_idx1 = np.squeeze(np.array(ligdents)).argsort()[-2:][::-1][0]
            bidentate_ligand_idx2 = np.squeeze(np.array(ligdents)).argsort()[-2:][::-1][1]
            bidentate_cons1 = ligcons[bidentate_ligand_idx1]
            bidentate_cons2 = ligcons[bidentate_ligand_idx2]
            pair_combos = list(combinations([0, 1, 2, 3], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in bidentate_cons1]+[mol.getAtom(ii).coords() for ii in bidentate_cons2])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            #### Seesaws will have only 1 ~180 degree angle, whereas planar ligands will have two. 
            #### Thus, after measuring the angles for all tetradentate connecting atoms through the metal,
            #### looking at the angle of the first (not zeroeth) element tells us whether or not we have a seesaw.
            test_angle = np.sort(np.array(angle_list))[::-1][1]
            if loud:
                print('ANGLE LIST',angle_list, 'sorted',np.sort(np.array(angle_list))[::-1])
            if test_angle < angle_cutoff:
                seesaw = True
                temp_cons = bidentate_cons1+bidentate_cons2
                axial_pair = [temp_cons[val] for val in list(pair_list[np.argmax(np.array(angle_list))])]
            else:
                seesaw = False
            if not seesaw:
                print('non seesaw')
                eq_lig_list = [bidentate_ligand_idx1,bidentate_ligand_idx2]
                eq_con_list = [bidentate_cons1,bidentate_cons2]  ### ADDED
                ax_lig_list = list(set(allowed)-set([bidentate_ligand_idx1,bidentate_ligand_idx2]))
                ax_con_list = [ligcons[j] for j in ax_lig_list] 
            else:  # 2 points in eq plane, seesaw case
                eq_ligcons = set([flat_ligcons[j] for j in eq_points])
                eq_con_bidentate_list = [list(set(ligcons[0]).intersection(eq_ligcons)),
                                         list(set(ligcons[1]).intersection(eq_ligcons)),
                                         list(set(ligcons[2]).intersection(eq_ligcons)),
                                         list(set(ligcons[3]).intersection(eq_ligcons))]
                eq_con_bidentate_list = [val for val in eq_con_bidentate_list if len(val)>0]
                eq_con_list = eq_con_bidentate_list
                flat_eq_con_list = [item for sublist in eq_con_list for item in sublist]
                eq_lig_list = [j for j, val in enumerate(ligcons) if len(set(val).intersection(set(flat_eq_con_list)))>0]
                ax_con_list = [[val] for val in flat_ligcons if val not in flat_eq_con_list]
                flat_ax_con_list = [item for sublist in ax_con_list for item in sublist]
                ax_lig_list = [j for j, val in enumerate(ligcons) if len(set(val).intersection(set(flat_ax_con_list)))>0]
    elif (n_ligs == 3):  # triple bidentate or 4+1+1, can be seesaw/planar or 3+2+1
        if max(ligdents) == 4:
            print('4+1+1 case')
            allowed = range(0, 3)
            tetradentate_ligand_idx = np.argmax(ligdents)
            tetradentate_cons = ligcons[tetradentate_ligand_idx]
            pair_combos = list(combinations([0, 1, 2, 3], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in tetradentate_cons])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            #### Seesaws will have only 1 ~180 degree angle, whereas planar ligands will have two. 
            #### Thus, after measuring the angles for all tetradentate connecting atoms through the metal,
            #### looking at the angle of the first (not zeroeth) element tells us whether or not we have a seesaw.
            test_angle = np.sort(np.array(angle_list))[::-1][1]
            if test_angle < angle_cutoff:
                seesaw = True
                #### In the seesaw, the two furthest apart are denoted as axial. 
                #### The equatorial plane consists of 2 seesaw connecting atoms, and two monodentates.
                axial_pair = [tetradentate_cons[val] for val in list(pair_list[np.argmax(np.array(angle_list))])]
            else:
                seesaw = False
                #### In the planar set, the two monodentates are axial, and tetradentate is equatorial.
                axial_pair = list(set(allowed) - set([tetradentate_ligand_idx]))
            if not seesaw:
                print('all tetradentate in eq plane')
                eq_lig_list = [tetradentate_ligand_idx]
                eq_con_list = [tetradentate_cons]  ### ADDED
                ax_lig_list = axial_pair
                ax_con_list = [ligcons[axial_pair[0]], ligcons[axial_pair[1]]]
            else:  # 2 points in eq plane, seesaw case
                print('two tetradentate in eq plane')
                eq_plane = list(set(flat_ligcons) - set(axial_pair))
                #### Find the two connecting atoms of the seesaw that lie in the equatorial plane
                tet_eq = list(set(eq_plane).intersection(set(tetradentate_cons)))
                eq_con_list = [ligcons[j] if len(ligcons[j]) == 1 else tet_eq for j in range(len(ligcons))]
                eq_lig_list = allowed
                ax_lig_list = [tetradentate_ligand_idx, tetradentate_ligand_idx]
                ax_con_list = [[axial_pair[0]], [axial_pair[1]]]
        elif max(ligdents) == 3:
            print('3+2+1 case')
            allowed = range(0, 3)
            tridentate_ligand_idx = np.argmax(ligdents)
            monodentate_ligand_idx = np.argmin(ligdents)
            bidentate_ligand_idx = list(set(allowed) - set([tridentate_ligand_idx]) - set([monodentate_ligand_idx]))[0]
            tridentate_cons = ligcons[tridentate_ligand_idx]
            bidentate_cons = ligcons[bidentate_ligand_idx]
            monodentate_cons = ligcons[monodentate_ligand_idx]
            pair_combos = list(combinations([0, 1, 2], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            #### The tridentate can either be planar or not. If planar, it will have at least one
            #### angle that is close to 180 between connecting atoms, through the metal.
            test_angle = np.sort(np.array(angle_list))[::-1][0]
            if test_angle < angle_cutoff:
                planar = False
                coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons])
                p1 = np.array(mol.getAtom(monodentate_cons[0]).coords())
                angle_list = []
                for k, tri_con in enumerate(tridentate_cons):
                    p2 = np.squeeze(np.array(coord_list[k]))
                    v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                    v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                    angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                    angle_list.append(angle)
                #### If not planar, one connection atom (for the monodentate ligand) will be opposite
                #### the connections for the tridentate ligands. Assigning this atom as the opposite monodentate.
                #### If not planar, considering the plane with the bidentate to be equatorial.
                opposite_monodentate = [tridentate_cons[np.argmax(np.array(angle_list))]]
            else:
                print('planar')
                planar = True
                coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons])
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                for bi_con in bidentate_cons:
                    p1 = np.array(mol.getAtom(bi_con).coords())
                    angle_list = []
                    for k, tri_con in enumerate(tridentate_cons):
                        p2 = np.squeeze(np.array(coord_list[k]))
                        v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                        v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                        angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                        angle_list.append(angle)
                    test_angle = np.sort(np.array(angle_list))[::-1][0]
                    #### If the tridentate is planar, only one of the bidentate connecting atoms
                    #### will be across the tridentate. Thus the equatorial plane will be defined
                    #### as the planar tridentate plus one connection atom for the bidentate
                    #### The monodentate and the second bidentate connecting atoms are denoted axial.
                    if test_angle > angle_cutoff:
                        bidentate_eq_con = [bi_con]
                        bidentate_ax_con = [list(set(bidentate_cons) - set([bi_con]))[0]]
            if planar:
                eq_lig_list = [tridentate_ligand_idx, bidentate_ligand_idx]
                eq_con_list = [tridentate_cons, bidentate_eq_con]  ### ADDED
                ax_lig_list = [monodentate_ligand_idx, bidentate_ligand_idx]
                ax_con_list = [monodentate_cons, bidentate_ax_con]
            else:
                eq_lig_list = [tridentate_ligand_idx, bidentate_ligand_idx]
                eq_con_list = [tridentate_cons, bidentate_cons]  ### ADDED
                ax_lig_list = [monodentate_ligand_idx, tridentate_ligand_idx]
                ax_con_list = [monodentate_cons, opposite_monodentate]
        elif (max(ligdents) == 2) and (min(ligdents) == 2):
            print('triple bidentate case')
            allowed = range(0, 3)
            bidentate_cons_1 = ligcons[0]
            bidentate_cons_2 = ligcons[1]
            bidentate_cons_3 = ligcons[2]
            #### Using previously defined eq plane to determine the axial positions
            eq_ligcons = set([flat_ligcons[j] for j in eq_points])
            eq_con_bidentate_list = [list(set(bidentate_cons_1).intersection(eq_ligcons)),
                                     list(set(bidentate_cons_2).intersection(eq_ligcons)),
                                     list(set(bidentate_cons_3).intersection(eq_ligcons))]
            #### all 3 ligands classified as equatorial because at least one atom connects to the equatorial plane for all
            #### axial only has two, where the axial connection is defined
            eq_con_list = [eq_con_bidentate_list[0], eq_con_bidentate_list[1], eq_con_bidentate_list[2]]
            eq_lig_list = allowed
            ax_con_list = [list(set(ligcons[i]) - set(eq_con_bidentate_list[i])) for i in allowed if
                           (len(eq_con_bidentate_list[i]) == 1)]
            ax_lig_list = [i for i in allowed if (len(eq_con_bidentate_list[i]) == 1)]
    elif (n_ligs == 2 and not pentadentate):  # 4+2 or 3+3
        if ((max(ligdents) == 4) and (min(ligdents) == 2)):
            print('4+2 mode')
            #### 4+2 cases are handled in the exact same way as 4+1+1 seesaws. See above for explanation.
            #### The seesaw tetradentate ligand has two axial and two equatorial connecting atoms.
            allowed = range(0, 2)
            tetradentate_ligand_idx = np.argmax(ligdents)
            tetradentate_cons = ligcons[tetradentate_ligand_idx]
            pair_combos = list(combinations([0, 1, 2, 3], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in tetradentate_cons])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            axial_pair = [tetradentate_cons[val] for val in list(pair_list[np.argmax(np.array(angle_list))])]
            eq_plane = list(set(flat_ligcons) - set(axial_pair))
            tet_eq = list(set(eq_plane).intersection(set(tetradentate_cons)))
            eq_con_list = [ligcons[j] if len(ligcons[j]) == 2 else tet_eq for j in range(len(ligcons))]
            eq_lig_list = allowed
            ax_lig_list = [tetradentate_ligand_idx, tetradentate_ligand_idx]
            ax_con_list = [[axial_pair[0]], [axial_pair[1]]]
        if ((max(ligdents) == 3) and (min(ligdents) == 3)):
            print('3+3 mode')
            allowed = range(0, 2)
            tridentate_cons_1 = ligcons[0]
            tridentate_cons_2 = ligcons[1]
            pair_combos = list(combinations([0, 1, 2], 2))
            angle_list = []
            pair_list = []
            coord_list = np.array([mol.getAtom(ii).coords() for ii in tridentate_cons_1])
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            test_angle = np.sort(np.array(angle_list))[::-1][0]

            #### Like the above 3+2+1 case, the tridentate must be classified into planar or not.
            #### If not planar, both ligands are both axial and equatorial.
            #### If planar, the ligand that falls on the eq_points (as found above), is only equatorial,
            #### but the tridentate ligand with one equatorial position and two axial is classified as both.
            if test_angle < angle_cutoff:
                planar = False
            else:
                planar = True
            eq_ligcons = set([flat_ligcons[j] for j in eq_points])
            if not planar:
                eq_con_list = [list(set(tridentate_cons_1).intersection(eq_ligcons)),
                               list(set(tridentate_cons_2).intersection(eq_ligcons))]
                ax_con_list = [list(set(tridentate_cons_1) - set(eq_con_list[0])),
                               list(set(tridentate_cons_2) - set(eq_con_list[1]))]
                eq_lig_list = [0, 1]
                ax_lig_list = [0, 1]
            else:
                eq_lig1 = list(set(tridentate_cons_1).intersection(eq_ligcons))
                eq_lig2 = list(set(tridentate_cons_2).intersection(eq_ligcons))
                eq_con_list = [eq_lig1, eq_lig2]
                eq_lig_list = [0, 1]
                if len(eq_lig1) == 1:
                    ax_con_list = [list(set(tridentate_cons_1) - set(eq_con_list[0]))]
                    ax_lig_list = [0]
                else:
                    ax_con_list = [list(set(tridentate_cons_2) - set(eq_con_list[1]))]
                    ax_lig_list = [1]
    elif (n_ligs == 2 and pentadentate):
        #### Handling for pentadentate scaffolds ####
        if loud:
            print('pentadentate case')
        allowed = [0, 1]
        not_eq = list()
        for j, built_ligs in enumerate(built_ligand_list):
            if len(ligcons[j]) == 1:
                #### This is the axial ligand ####
                # print(j, 'axial lig')
                top_lig = j
                top_con = ligcons[j]
                not_eq.append(top_lig)
            else:
                pentadentate_coord_list = np.array([mol.getAtom(
                    ii).coords() for ii in ligcons[j]])
                ##### Adjusting this so that by default, any 4 within the same plane will be assigned as eq. ###
                if loud:
                    print('pentadentate coord LIST!')
                    print(pentadentate_coord_list)
                point_combos = combinations([0, 1, 2, 3, 4], 4)
                error_list = []
                combo_list = []
                for i, combo in enumerate(point_combos):
                    combo_list.append(list(combo))
                    A = []
                    b = []
                    for point_num in combo:
                        coordlist = pentadentate_coord_list[point_num]
                        A.append([coordlist[0], coordlist[1], 1])
                        b.append(coordlist[2])
                    ##### This code builds the best fit plane between 4 points,
                    ##### Then calculates the variance of the 4 points with respect to the plane
                    ##### The 4 that have the least variance are flagged as the eq plane.
                    mat_b = np.matrix(b).T
                    mat_A = np.matrix(A)
                    try:
                        fit = (mat_A.T * mat_A).I * mat_A.T * mat_b
                        errors = np.squeeze(np.array(mat_b - mat_A * fit))
                        error_var = np.var(errors)
                        error_list.append(error_var)
                    except:
                        #### This will catch if it errors due to a singular matrix. Means the fit is perfect.
                        #### If the fit is perfect, then it is the equatorial plane. Assigned this way.
                        error_list.append(0)
                if loud:
                    print('combos below')
                    print(combo_list)
                    print('errors next, argmin combo selected')
                    print(error_list)
                not_ax_points = combo_list[np.argmin(error_list)]
                if len(set(not_ax_points)) != 4:
                    print('The equatorial plane is not being assigned correctly. Please check.')
                    sardines
                else:
                    bot_idx = list(set(range(5)) - set(not_ax_points))[0]
                    if loud:
                        print('This is bot_idx', bot_idx)
                    bot_lig = j
                    bot_con = [ligcons[j][bot_idx]]

        allowed = list(set(allowed) - set(not_eq))
        if loud:
            print('this is the allowed list', allowed, not_eq)
        eq_lig_list = allowed
        eq_con_list = [
            list(set([ligcons[i] for i in allowed][0]) - set(top_con) - set(bot_con))]
        ax_lig_list = [top_lig, bot_lig]
        ax_con_list = [top_con, bot_con]
        if loud:
            print('con lists', eq_con_list, ax_con_list)
        ###########################################################################################
        # In the above, the pentadentate ligand is classified as both axial and equatorial.       #
        # The lc atoms are decided by the z-position. Thus the pentadentate ligand has 4 eq-lc    #
        # and 1 ax-lc. Currently should be able to check this and set that up.                    #
        ###########################################################################################
    elif (n_ligs == 1 and hexadentate):
        allowed = [0, 1]
        if loud:
            print('hexadentate case')
        not_eq = list()
        for j, built_ligs in enumerate(built_ligand_list):
            hexadentate_coord_list = np.array([mol.getAtom(
                ii).coords() for ii in ligcons[j]])
            ##### Adjusting this so that by default, any 4 within the same plane will be assigned as eq. ###
            if loud:
                print('hexadentate coord LIST!')
                print(hexadentate_coord_list)
            # point_combos = combinations([1,2,3,4,5,6],4)
            pair_combos = list(combinations([0, 1, 2, 3, 4, 5], 2))
            angle_list = []
            pair_list = []
            for i, pair in enumerate(pair_combos):
                pair_list.append(list(pair))
                p1 = np.squeeze(np.array(hexadentate_coord_list[list(pair)[0]]))
                p2 = np.squeeze(np.array(hexadentate_coord_list[list(pair)[1]]))
                m = np.array([mol.getAtom(mol.findMetal()[0]).coords()])
                v1u = np.squeeze(np.array((m - p1) / np.linalg.norm((m - p1))))
                v2u = np.squeeze(np.array((m - p2) / np.linalg.norm((m - p2))))
                # print('v1v2',v1u,v2u)
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))
                if loud:
                    print('pair of atoms, then angle', pair, angle)
                angle_list.append(angle)
            argsort_angle_list = np.squeeze(np.array(angle_list)).argsort()[-3:][::-1]
            point_combos = [pair_list[argsort_angle_list[0]] + pair_list[argsort_angle_list[1]],
                            pair_list[argsort_angle_list[1]] + pair_list[argsort_angle_list[2]],
                            pair_list[argsort_angle_list[2]] + pair_list[argsort_angle_list[0]]]
            error_list = []
            combo_list = []
            fitlist = []
            for i, combo in enumerate(point_combos):
                combo_list.append(combo)
                A = []
                b = []
                for point_num in combo:
                    coordlist = hexadentate_coord_list[point_num]
                    A.append([coordlist[0], coordlist[1], 1])
                    b.append(coordlist[2])
                ##### This code builds the best fit plane between 4 points,
                ##### Then calculates the variance of the 4 points with respect to the plane
                ##### The 4 that have the least variance are flagged as the eq plane.
                mat_b = np.matrix(b).T
                mat_A = np.matrix(A)
                try:
                    fit = (mat_A.T * mat_A).I * mat_A.T * mat_b
                    errors = np.squeeze(np.array(mat_b - mat_A * fit))
                    error_var = np.var(errors)
                    error_list.append(error_var)
                except:
                    #### This will catch if it errors due to a singular matrix. Means the fit is perfect.
                    #### If the fit is perfect, then it is the equatorial plane. Assigned this way.
                    error_list.append(0)
            if loud:
                print('combos below')
                print(combo_list)
                print('errors next, argmin combo selected')
                print(error_list)
            best_fit_planes = np.squeeze(np.array(error_list)).argsort()[:3]
            perpdist = []
            perpcombo = []
            for fitnum, best_fit in enumerate(best_fit_planes):
                temp_fit = fitlist[best_fit]
                temp_combo = combo_list[best_fit]
                perpcombo.append(int(best_fit))
                temp_ax = set(range(0, 6)) - set(temp_combo)
                ax_dist = []
                for point_num in temp_ax:
                    coordlist = hexadentate_coord_list[point_num]
                    planez = [coordlist[0], coordlist[1], 1] * temp_fit
                    # planez = temp_fit[0] * coordlist[0] + temp_fit[1] * coordlist[1] + fit[2]
                    plane_coords = [coordlist[0], coordlist[1], planez]
                    adjusted_coords = [coordlist[0], coordlist[1], coordlist[2]]
                    squared_dist = np.sum((np.array(adjusted_coords) - np.array(plane_coords)) ** 2)
                    dist = np.squeeze(np.sqrt(squared_dist))
                    if loud:
                        print('dist', dist)
                    ax_dist.append(dist)
                perpdist.append(np.mean(ax_dist))
            if loud:
                print("Perpendicular distance is", perpdist, perpcombo, len(perpdist), len(best_fit_planes))
            not_ax_points = combo_list[perpcombo[np.argmax(np.array(perpdist))]]
            if len(set(not_ax_points)) != 4:
                print('The equatorial plane is not being assigned correctly. Please check.')
                sardines
            else:
                bot_idx = list(set(range(6)) - set(not_ax_points))[0]
                top_idx = list(set(range(6)) - set(not_ax_points))[1]
                if loud:
                    print('This is bot_idx', bot_idx)
                bot_lig = j
                top_lig = j
                bot_con = [ligcons[j][bot_idx]]
                top_con = [ligcons[j][top_idx]]
        allowed = list(set(allowed) - set(not_eq))
        if loud:
            print('this is the allowed list', allowed, not_eq)
        eq_lig_list = [top_lig]
        eq_con_list = [
            list(set([ligcons[0][i] for i in not_ax_points]))]
        ax_lig_list = [top_lig, bot_lig]
        ax_con_list = [top_con, bot_con]
        if loud:
            print('con lists', eq_con_list, ax_con_list)

    ax_ligand_list = [built_ligand_list[i] for i in ax_lig_list]
    eq_ligand_list = [built_ligand_list[i] for i in eq_lig_list]
    if loud and valid:
        print(('lig_nat_list', lig_natoms_list))
        print(('eq_liq is ind ', eq_lig_list))
        print(('ax_liq is ind ', ax_lig_list))
        print(('ax built lig [0] ext ind :' +
               str(list(built_ligand_list[ax_lig_list[0]].ext_int_dict.keys()))))
        if len(ax_lig_list) > 1:
            print(('ax built lig [1] ext ind :' +
                   str(list(built_ligand_list[ax_lig_list[1]].ext_int_dict.keys()))))
        print(('eq built lig [0] ext ind: ' +
               str(list(built_ligand_list[eq_lig_list[0]].ext_int_dict.keys()))))
        print(('eq_con is ' + str((eq_con_list))))
        print(('ax_con is ' + str((ax_con_list))))
    if name:  ## TODO: this part might be broken. Maybe not of much use though.
        for i, ax_ligand in enumerate(ax_ligand_list):
            if not os.path.isdir('ligands'):
                os.mkdir('ligands')
            ax_ligand.mol.writexyz('ligands/' + name +
                                   '_' + str(i) + '_ax.xyz')
        for i, eq_ligand in enumerate(eq_ligand_list):
            if not os.path.isdir('ligands'):
                os.mkdir('ligands')
            eq_ligand.mol.writexyz('ligands/' + name +
                                   '_' + str(i) + '_eq.xyz')
    for j, ax_con in enumerate(ax_con_list):
        current_ligand_index_list = built_ligand_list[ax_lig_list[j]].index_list
        ax_con_int_list.append([current_ligand_index_list.index(i) for i in ax_con])
        # ax_con_int_list.append(
        #     [built_ligand_list[ax_lig_list[j]].ext_int_dict[i] for i in ax_con])  # convert to interal index ## old one with obtain_mol3d
    for j, eq_con in enumerate(eq_con_list):
        current_ligand_index_list = built_ligand_list[eq_lig_list[j]].index_list
        eq_con_int_list.append([current_ligand_index_list.index(i) for i in eq_con])
        # eq_con_int_list.append(
        #     [built_ligand_list[eq_lig_list[j]].ext_int_dict[i] for i in eq_con])  # convert to interal index ## old one with obtain_mol3d
    if loud:
        print(('int eq ' + str(eq_con_int_list)))
        print(('ext eq ' + str(eq_con_list)))
        print('**********************************************')
    for ax_lig in ax_lig_list:
        ax_natoms_list.append(lig_natoms_list[ax_lig])
    for eq_lig in eq_lig_list:
        eq_natoms_list.append(lig_natoms_list[eq_lig])
    return ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, built_ligand_list
