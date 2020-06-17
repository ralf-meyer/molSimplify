# Written JP Janet
# for HJK Group
# Dpt of Chemical Engineering, MIT

##########################################################
######## Defines methods for assembling    ###############
########     RACs from lists of ligands    ###############
##########################################################
from __future__ import print_function
import numpy as np
import sys
from molSimplify.Classes.ligand import ligand_breakdown
from molSimplify.Classes.ligand import ligand_assign_consistent as ligand_assign
from molSimplify.Classes.globalvars import globalvars

globs = globalvars()

## Gets the RACs of an octahedral complex without/without geo
#  @param this_complex mol3D() we want RACs for
#  @param custom_ligand_dict optional dict defining ligands (see below)
#  @return descriptor_names updated names
#  @return descriptors updated RACs
def get_descriptor_vector(this_complex,custom_ligand_dict=False,ox_modifier=False, NumB=False, Zeff=False, \
    lacRACs = True, loud = False):
    descriptor_names = []
    descriptors = []
    # Generate custom_ligand_dict if one not passed!
    if not custom_ligand_dict:
        if lacRACs:
            from molSimplify.Classes.ligand import ligand_assign_consistent as ligand_assign
        else:
            from molSimplify.Classes.ligand import ligand_assign as ligand_assign
        liglist, ligdents, ligcons = ligand_breakdown(this_complex)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, \
            ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, \
                built_ligand_list = ligand_assign(
                this_complex, liglist, ligdents, ligcons, loud, name=False)
        custom_ligand_dict = {'ax_ligand_list':ax_ligand_list, 'eq_ligand_list':eq_ligand_list,
                               'ax_con_int_list':ax_con_int_list, 'eq_con_int_list':eq_con_int_list}
    ## misc descriptors
    results_dictionary = generate_all_ligand_misc(this_complex,loud=False,
                                                    custom_ligand_dict=custom_ligand_dict)
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_ax'],'misc','ax')
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_eq'],'misc','eq')

    ## full ACs
    results_dictionary = generate_full_complex_autocorrelations(this_complex,depth=3,loud=False,flag_name=False,
                                                                modifier=ox_modifier, NumB=NumB, Zeff=Zeff)
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['results'],'f','all')

    ## ligand ACs
    #print('get ligand ACs')
    results_dictionary = generate_all_ligand_autocorrelations(this_complex,depth=3,loud=False,name=False,
                                                                flag_name=False,
                                                                custom_ligand_dict=custom_ligand_dict,
                                                                NumB=NumB, Zeff=Zeff)
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_ax_full'],'f','ax')
    descriptor_names, descriptors =  append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_eq_full'],'f','eq')
    descriptor_names, descriptors =  append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_ax_con'],'lc','ax')
    descriptor_names, descriptors =  append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_eq_con'],'lc','eq')

    results_dictionary = generate_all_ligand_deltametrics(this_complex,depth=3,loud=False,name=False,
                                                            custom_ligand_dict=custom_ligand_dict,
                                                            NumB=NumB, Zeff=Zeff)
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_ax_con'],'D_lc','ax')
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['result_eq_con'],'D_lc','eq')

    ## metal ACs
    #print('getting metal ACs')
    results_dictionary = generate_metal_autocorrelations(this_complex,depth=3,loud=False,
                                                            modifier=ox_modifier,
                                                            NumB=NumB,Zeff=Zeff)
    descriptor_names, descriptors =  append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['results'],'mc','all')

    results_dictionary = generate_metal_deltametrics(this_complex,depth=3,loud=False,
                                                        modifier=ox_modifier,
                                                        NumB=NumB,Zeff=Zeff)
    descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['results'],'D_mc','all')

    # ## ox-metal ACs, if ox available
    if ox_modifier:
        results_dictionary = generate_metal_ox_autocorrelations(ox_modifier, this_complex,depth=3,loud=False)
        descriptor_names, descriptors =  append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['results'],'mc','all')
        results_dictionary = generate_metal_ox_deltametrics(ox_modifier,this_complex,depth=3,loud=False)
        descriptor_names, descriptors = append_descriptors(descriptor_names, descriptors,
                                                        results_dictionary['colnames'],results_dictionary['results'],'D_mc','all')
    return descriptor_names, descriptors


## Gets the derivatives of RACs of an octahedral complex with geo only!
#  @param this_complex mol3D() we want derivatives for
#  @param custom_ligand_dict optional dict defining ligands (see below)
#  @return descriptor_derivative_names matrix of names
#  @return descriptor_derivatives derivatives of racs w.r.t atomic props
def get_descriptor_derivatives(this_complex, custom_ligand_dict=False, ox_modifier=False, lacRACs=True, loud=False):
    if not custom_ligand_dict:
        if lacRACs:
            from molSimplify.Classes.ligand import ligand_assign_consistent as ligand_assign
        else:
            from molSimplify.Classes.ligand import ligand_assign as ligand_assign
        liglist, ligdents, ligcons = ligand_breakdown(this_complex)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, \
            ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, \
                built_ligand_list = ligand_assign(
                this_complex, liglist, ligdents, ligcons, loud, name=False)
        custom_ligand_dict = {'ax_ligand_list':ax_ligand_list, 'eq_ligand_list':eq_ligand_list,
                               'ax_con_int_list':ax_con_int_list, 'eq_con_int_list':eq_con_int_list}
    ##  cannot do misc descriptors !
    descriptor_derivative_names = []
    descriptor_derivatives = None
    ## full ACs
    results_dictionary = generate_full_complex_autocorrelation_derivatives(this_complex,depth=3,
                                                                           loud=False,flag_name=False,
                                                                           modifier=ox_modifier)
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['results'],'f','all')
    ## ligand ACs
    #print('getting ligand AC derivatives')
    results_dictionary = generate_all_ligand_autocorrelation_derivatives(this_complex,depth=3,loud=False,name=False, custom_ligand_dict=custom_ligand_dict)
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_ax_full'],'f','ax')
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_eq_full'],'f','eq')
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_ax_con'],'lc','ax')
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_eq_con'],'lc','eq')
    results_dictionary = generate_all_ligand_deltametric_derivatives(this_complex,depth=3,loud=False,name=False, custom_ligand_dict=custom_ligand_dict)
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_ax_con'],'D_lc','ax')
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['result_eq_con'],'D_lc','eq')
    ## metal ACs
    #print('getting metal AC derivatives')
    results_dictionary = generate_metal_autocorrelation_derivatives(this_complex,depth=3,loud=False,modifier=ox_modifier)
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['results'],'mc','all')
    results_dictionary = generate_metal_deltametric_derivatives(this_complex,depth=3,loud=False,modifier=ox_modifier)
    descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                        results_dictionary['colnames'],results_dictionary['results'],'D_mc','all')
    # ## ox-metal ACs
    if ox_modifier:
        results_dictionary = generate_metal_ox_autocorrelation_derivatives(ox_modifier, this_complex,depth=3,loud=False)
        descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                            results_dictionary['colnames'],results_dictionary['results'],'mc','all')

        results_dictionary =  generate_metal_ox_deltametric_derivatives(ox_modifier, this_complex,depth=3,loud=False)
        descriptor_derivative_names, descriptor_derivatives = append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives,
                                                                                            results_dictionary['colnames'],results_dictionary['results'],'D_mc','all')

    return descriptor_derivative_names, descriptor_derivatives


## utility to build standardly formated RACS
#  @param descriptor_names RAC names, will be appended to
#  @param descriptors RAC, will be appended to
#  @param list_of_names names, will be added
#  @param list_of_props types of RACs
#  @param prefix RAC prefix
#  @param suffix RAC suffix
#  @return descriptor_names updated names
#  @return descriptors updated RACs
def append_descriptors(descriptor_names,descriptors,list_of_names,list_of_props,prefix,suffix):
    try:
        basestring
    except NameError:
        basestring = str

    for names in list_of_names:
        if not isinstance(names, basestring):
            names = ["-".join([prefix,str(i),suffix]) for i in names]
            descriptor_names += names
        else:
            names = "-".join([prefix,str(names),suffix])
            descriptor_names.append(names)
    for values in list_of_props:
        if not isinstance(names, basestring):
            descriptors.extend(values)
        else:
            descriptors.append(values)
    return descriptor_names, descriptors

## utility to build standardly formated RACS derivatives
#  @param descriptor_derivative_names RAC names, will be a matrix!
#  @param descriptor_derivatives RAC, will be appended to
#  @param mat_of_names names, will be added
#  @param dmat mat of RACs
#  @param prefix RAC prefix
#  @param suffix RAC suffix
#  @return descriptor_derivative_names updated names
#  @return descriptor_derivatives updated RACs
def append_descriptor_derivatives(descriptor_derivative_names,descriptor_derivatives, mat_of_names,dmat,prefix,suffix):
    try:
        basestring
    except NameError:
        basestring = str

    for names in mat_of_names:
        jnames = ["-".join([prefix,str(i),suffix]) for i in names]
        descriptor_derivative_names.append(jnames)
    if descriptor_derivatives is None:
        descriptor_derivatives = dmat
    else:
        descriptor_derivatives =  np.row_stack([descriptor_derivatives,dmat])
    return descriptor_derivative_names, descriptor_derivatives

## this function returns the products autocorrelation 
## for one atom
#  @ param mol,  mol3D class
#  @ param prop_vec, vector, property of atoms in mol in order of index
#  @ param orig, int, zero-indexed starting atom
#  @ param d, int, number of hops to travel
#  @ param oct, bool, if complex is octahedral, will use better bond checks
#  @ return results_vector, vector, assembled autocorrelation
def autocorrelation(mol, prop_vec, orig, d, oct=True, catoms=None, use_dist=False):
    result_vector = np.zeros(d + 1)
    hopped = 0
    active_set = set([orig])
    historical_set = set()
    if not use_dist:
        result_vector[hopped] = prop_vec[orig] * prop_vec[orig]
    else:
        result_vector[hopped] = 0.5 * abs(prop_vec[orig]) ** 2.4 / mol.natoms
    while hopped < (d):
        hopped += 1
        new_active_set = set()
        for this_atom in active_set:
            this_atoms_neighbors = mol.getBondedAtomsSmart(this_atom, oct=oct)
            for bound_atoms in this_atoms_neighbors:
                if (bound_atoms not in historical_set) and (bound_atoms not in active_set):
                    new_active_set.add(bound_atoms)
        for inds in new_active_set:
            if not use_dist:
                result_vector[hopped] += prop_vec[orig] * prop_vec[inds]
            else:
                this_dist = mol.getDistToMetal(orig,inds)
                result_vector[hopped] += prop_vec[orig] * prop_vec[inds] / (this_dist * mol.natoms)
            historical_set.update(active_set)
        active_set = new_active_set
    return (result_vector)


## this function returns the derivative vector
## of the scalar autocorrelation 
## starting at orig with depth d,
## with respect to the atomic properties
## in prop_vec, for all atoms.
## The return type is np.array
## Be sure to read this carefully!
## for one atom
#  @ param mol,  mol3D class
#  @ param prop_vec, vector, property of atoms in mol in order of index
#  @ param orig, int, zero-indexed starting atom
#  @ param d, int, number of hops to travel
#  @ param oct, bool, if complex is octahedral, will use better bond checks
#  @ return derivative_mat, vector, autocorrelation_derivatives
def autocorrelation_derivative(mol, prop_vec, orig, d, oct=True, catoms=None):
    derivative_mat = np.zeros((d + 1, len(prop_vec)))
    # loop for each atom
    hopped = 0
    active_set = set([orig])
    historical_set = set()
    for derivate_ind in range(0, len(prop_vec)):
        if derivate_ind == orig:
            derivative_mat[hopped, derivate_ind] = 2 * prop_vec[orig]
        else:
            derivative_mat[hopped, derivate_ind] = 0
    while hopped < (d):

        hopped += 1
        new_active_set = set()
        for this_atom in active_set:
            ## prepare all atoms attached to this connection
            # print('called in AC')
            this_atoms_neighbors = mol.getBondedAtomsSmart(this_atom, oct=oct)
            for bound_atoms in this_atoms_neighbors:
                if (bound_atoms not in historical_set) and (bound_atoms not in active_set):
                    new_active_set.add(bound_atoms)
        # print('new active set at hop = ' +str(hopped) + ' is ' +str(new_active_set))
        for inds in new_active_set:
            for derivate_ind in range(0, len(prop_vec)):
                if derivate_ind == orig:
                    derivative_mat[hopped, derivate_ind] += prop_vec[inds]
                elif derivate_ind == inds:
                    derivative_mat[hopped, derivate_ind] += prop_vec[orig]
            historical_set.update(active_set)
        active_set = new_active_set
    return (derivative_mat)


## this function returns the deltametric autocorrelation
## over the whole molecule
# Inputs:
# @param mol, mol3D class
# @param prop_vec, vector, property of atoms in mol in order of index
# @param orig, int, zero-indexed starting atom
# @param d, int, number of hops to travel
# @param oct, bool, if complex is octahedral, will use better bond checks
# @return results_vector, vector, deltametric autocorrelations.
def deltametric(mol, prop_vec, orig, d, oct=True, catoms=None):
    result_vector = np.zeros(d + 1)
    hopped = 0
    active_set = set([orig])
    historical_set = set()
    result_vector[hopped] = 0.00
    while hopped < (d):
        hopped += 1
        new_active_set = set()
        for this_atom in active_set:
            ## prepare all atoms attached to this connection
            # print('called in DAC')
            this_atoms_neighbors = mol.getBondedAtomsSmart(this_atom, oct=oct)
            for bound_atoms in this_atoms_neighbors:
                if (bound_atoms not in historical_set) and (bound_atoms not in active_set):
                    new_active_set.add(bound_atoms)
        # print('new active set at hop = ' +str(hopped) + ' is ' +str(new_active_set))
        for inds in new_active_set:
            result_vector[hopped] += prop_vec[orig] - prop_vec[inds]
            historical_set.update(active_set)
        active_set = new_active_set
    return (result_vector)


## this function returns the derivative vector
## of the scalar deltametric autocorrelation
## starting at orig with depth d,
## with respect to the atomic properties
## in prop_vec, for all atoms.
## The return type is np.array for
## Be sure to read this carefully!
#  @ param mol,  mol3D class
#  @ param prop_vec, vector, property of atoms in mol in order of index
#  @ param orig, int, zero-indexed starting atom
#  @ param d, int, number of hops to travel
#  @ param oct, bool, if complex is octahedral, will use better bond checks
#  @ return derivative_mat, vector, autocorrelation_derivatives
def deltametric_derivative(mol, prop_vec, orig, d, oct=True, catoms=None):
    derivative_mat = np.zeros((d + 1, len(prop_vec)))
    hopped = 0
    active_set = set([orig])
    historical_set = set()
    ## the zero-depth element is always zero
    for derivate_ind in range(0, len(prop_vec)):
        derivative_mat[hopped, derivate_ind] = 0.0
    while hopped < (d):
        hopped += 1
        new_active_set = set()
        for this_atom in active_set:
            ## prepare all atoms attached to this connection
            # print('called in DAC')
            this_atoms_neighbors = mol.getBondedAtomsSmart(this_atom, oct=oct)
            for bound_atoms in this_atoms_neighbors:
                if (bound_atoms not in historical_set) and (bound_atoms not in active_set):
                    new_active_set.add(bound_atoms)
        # print('new active set at hop = ' +str(hopped) + ' is ' +str(new_active_set))
        for inds in new_active_set:
            for derivate_ind in range(0, len(prop_vec)):
                if derivate_ind == orig:
                    derivative_mat[hopped, derivate_ind] += 1
                elif derivate_ind == inds:
                    derivative_mat[hopped, derivate_ind] += -1
        historical_set.update(active_set)
        active_set = new_active_set
    return (derivative_mat)

## assigns the value of property
## for atom i (zero index) in mol
## to position i in returned vector
## can be used to create weighted
## graph representations
## oct - bool, if complex is octahedral, will use better bond checks
## modifier - dict, used to modify prop vector (e.g. for adding
##             ONLY used with  ox_nuclear_charge    ox or charge)
##              {"Fe":2, "Co": 3} etc
## note that ident just codes every atom as one, this gives
## a purely toplogical index. coord gives the number of
## connecting atom to attom i (similar to Randic index)
# Acceptable prop valuse: ['electronegativity', 'nuclear_charge', 'ident', 'topology',
# 'ox_nuclear_charge', 'size', 'vdwrad', 'effective_nuclear_charge', 'polarizability',
# 'bondvalence', 'num_bonds', 'bondvalence_devi', 'bodavrg', 'bodstd', 'charge']
# @param mol, mol3D class
# @param prop, str, property to assemble 
# @return w, vector, property vector for mol by atom.
def construct_property_vector(mol, prop, oct=True, modifier=False):
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology',
                       'ox_nuclear_charge', 'size', 'vdwrad', 'effective_nuclear_charge', 'polarizability',
                       'bondvalence', 'num_bonds', 'bondvalence_devi', 'bodavrg', 'bodstd', 'charge']
    prop_dict = dict()
    w = np.zeros(mol.natoms)
    done = False
    if not prop in allowed_strings:
        print(('error, property  ' + str(prop) + ' is not a vaild choice'))
        print((' options are  ' + str(allowed_strings)))
        return False
    if prop == 'electronegativity':
        prop_dict = globs.endict()
    elif prop == 'size':
        at_keys = list(globs.amass().keys())
        for keys in at_keys:
            values = globs.amass()[keys][2]
            prop_dict.update({keys: values})
    elif prop == 'nuclear_charge':
        at_keys = list(globs.amass().keys())
        for keys in at_keys:
            values = globs.amass()[keys][1]
            prop_dict.update({keys: values})
    elif prop == 'effective_nuclear_charge':  # Uses number of valence electrons
        at_keys = list(globs.amass().keys())
        for keys in at_keys:
            values = globs.amass()[keys][0]
            prop_dict.update({keys: values}) 
    elif prop == 'ox_nuclear_charge':
        if not modifier:
            print('Error, must give modifier with ox_nuclear_charge')
            return False
        else:
            at_keys = list(globs.amass().keys())
            for keys in at_keys:
                values = globs.amass()[keys][1]
                if keys in list(modifier.keys()):
                    values -= float(modifier[keys]) # assumes oxidation state provided (i.e. Fe(IV))
                prop_dict.update({keys: values})
    elif prop == 'polarizability':
        prop_dict =  globs.polarizability()
        for i, atoms in enumerate(mol.getAtoms()):
            atom_type = atoms.symbol()
            w[i] = prop_dict[atom_type]
    elif prop == 'ident':
        at_keys = list(globs.amass().keys())
        for keys in at_keys:
            prop_dict.update({keys: 1})
    elif prop == 'topology':
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = len(mol.getBondedAtomsSmart(i, oct=oct))
        done = True
    elif prop == 'vdwrad':
        prop_dict = globs.vdwrad()
        for i, atoms in enumerate(mol.getAtoms()):
            atom_type = atoms.symbol()
            if atom_type in globs.metalslist():
                w[i] = globs.amass()[atoms.symbol()][2]
            else:
                w[i] = prop_dict[atoms.symbol()]
        done = True
    elif prop == 'bondvalence':
        assert len(mol.getAtoms()) == len(mol.bv_dict)
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = mol.bv_dict[i]
        done = True
    elif prop == 'num_bonds':
        for i, atom in enumerate(mol.getAtoms()):
            if not atom.ismetal():
                w[i] = globs.bondsdict()[atom.symbol()]
            else:
                w[i] = len(mol.getBondedAtomsSmart(i, oct=True))
        done = True
    elif prop == 'bondvalence_devi':
        assert len(mol.getAtoms()) == len(mol.bvd_dict)
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = mol.bvd_dict[i]
        done = True
    elif prop == 'bodavrg':
        assert len(mol.getAtoms()) == len(mol.bodavrg_dict)
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = mol.bodavrg_dict[i]
        done = True
    elif prop == 'bodstd':
        assert len(mol.getAtoms()) == len(mol.bodstd_dict)
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = mol.bodstd_dict[i]
        done = True
    elif prop == 'charge':
        assert len(mol.getAtoms()) == len(mol.charge_dict)
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = mol.charge_dict[i]
        done = True
    if not done:
        for i, atoms in enumerate(mol.getAtoms()):
            w[i] = prop_dict[atoms.symbol()]
    return (w)


def full_autocorrelation(mol, prop, d, oct=oct, modifier=False, use_dist=False):
    w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
    index_set = list(range(0, mol.natoms))
    autocorrelation_vector = np.zeros(d + 1)
    for centers in index_set:
        autocorrelation_vector += autocorrelation(mol, w, centers, d, oct=oct, use_dist=use_dist)
    return (autocorrelation_vector)


def full_autocorrelation_derivative(mol, prop, d, oct=oct, modifier=False):
    w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
    index_set = list(range(0, mol.natoms))
    autocorrelation_derivative_mat = np.zeros((d + 1, mol.natoms))
    for centers in index_set:
        autocorrelation_derivative_mat += autocorrelation_derivative(mol, w, centers, d, oct=oct)
    return (autocorrelation_derivative_mat)


def generate_full_complex_autocorrelations(mol, loud,
                                           depth=4, oct=True,
                                           flag_name=False, modifier=False,
                                           use_dist=False, NumB=False, Zeff=False, polarizability=False):
    result = list()
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    if polarizability:
        allowed_strings += ["polarizability"]
        labels_strings += ["alpha"]
    for ii, properties in enumerate(allowed_strings):
        metal_ac = full_autocorrelation(mol, properties, depth,
                                        oct=oct, modifier=modifier,
                                        use_dist=use_dist)
        this_colnames = []
        for i in range(0, depth + 1):
            this_colnames.append(labels_strings[ii] + '-' + str(i))
        colnames.append(this_colnames)
        result.append(metal_ac)
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_f_all': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary

def generate_full_complex_autocorrelation_derivatives(mol, loud, depth=4, oct=True, flag_name=False,
                                                      modifier=False, NumB=False, Zeff=False):
    result = None
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    for ii, properties in enumerate(allowed_strings):
        f_ac_der = full_autocorrelation_derivative(mol, properties, depth, oct=oct, modifier=modifier)
        for i in range(0, depth + 1):
            colnames.append(['d' + labels_strings[ii] + '-' + str(i) + '/d' + labels_strings[ii] + str(j) for j in
                             range(0, mol.natoms)])
        # colnames.append(this_colnames)
        if result is None:
            result = f_ac_der
        else:
            result = np.row_stack([result, f_ac_der])
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_f_all': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def atom_only_autocorrelation(mol, prop, d, atomIdx, oct=True):
    ## atomIdx must b either a list of indcies
    ## or a single index
    w = construct_property_vector(mol, prop, oct)
    autocorrelation_vector = np.zeros(d + 1)
    if hasattr(atomIdx, "__len__"):
        for elements in atomIdx:
            autocorrelation_vector += autocorrelation(mol, w, elements, d, oct=oct)
        autocorrelation_vector = np.divide(autocorrelation_vector, len(atomIdx))
    else:
        autocorrelation_vector += autocorrelation(mol, w, atomIdx, d, oct=oct)
    return (autocorrelation_vector)


def atom_only_autocorrelation_derivative(mol, prop, d, atomIdx, oct=True):
    ## atomIdx must b either a list of indcies
    ## or a single index
    w = construct_property_vector(mol, prop, oct)
    autocorrelation_derivative_mat = np.zeros((d + 1, mol.natoms))
    if hasattr(atomIdx, "__len__"):
        for elements in atomIdx:
            autocorrelation_derivative_mat += autocorrelation_derivative(mol, w, elements, d, oct=oct)
        autocorrelation_derivative_mat = np.divide(autocorrelation_derivative_mat, len(atomIdx))
    else:
        autocorrelation_derivative_mat += autocorrelation_derivative(mol, w, atomIdx, d, oct=oct)
    return (autocorrelation_derivative_mat)


def metal_only_autocorrelation(mol, prop, d, oct=True, catoms=None,
                               func=autocorrelation, modifier=False):
    autocorrelation_vector = np.zeros(d)
    try:
        metal_ind = mol.findMetal()[0]
        w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
        autocorrelation_vector = func(mol, w, metal_ind, d, oct=oct,
                                      catoms=catoms)
    except:
        print('Error, no metal found in mol object!')
        return False
    return (autocorrelation_vector)


def metal_only_autocorrelation_derivative(mol, prop, d, oct=True, catoms=None,
                                          func=autocorrelation_derivative, modifier=False):
    autocorrelation_vector_derivative = np.zeros((d + 1, len(prop)))
    try:
        metal_ind = mol.findMetal()[0]
        w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
        autocorrelation_vector_derivative = func(mol, w, metal_ind, d, oct=oct,
                                                 catoms=catoms)
    except:
        print('Error, no metal found in mol object!')
        return False
    return (autocorrelation_vector_derivative)


def atom_only_deltametric(mol, prop, d, atomIdx, oct=True, modifier=False):
    ## atomIdx must b either a list of indcies
    ## or a single index
    w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)

    deltametric_vector = np.zeros(d + 1)
    if hasattr(atomIdx, "__len__"):
        for elements in atomIdx:
            deltametric_vector += deltametric(mol, w, elements, d, oct=oct)
        deltametric_vector = np.divide(deltametric_vector, len(atomIdx))
    else:
        deltametric_vector += deltametric(mol, w, atomIdx, d, oct=oct)
    return (deltametric_vector)


def atom_only_deltametric_derivative(mol, prop, d, atomIdx, oct=True, modifier=False):
    ## atomIdx must b either a list of indcies
    ## or a single index
    w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)

    deltametric_derivative_mat = np.zeros((d + 1, mol.natoms))
    if hasattr(atomIdx, "__len__"):
        for elements in atomIdx:
            deltametric_derivative_mat += deltametric_derivative(mol, w, elements, d, oct=oct)
        deltametric_derivative_mat = np.divide(deltametric_derivative_mat, len(atomIdx))
    else:

        deltametric_derivative_mat += deltametric_derivative(mol, w, atomIdx, d, oct=oct)
    return (deltametric_derivative_mat)


def metal_only_deltametric_derivative(mol, prop, d, oct=True, catoms=None,
                                      func=deltametric_derivative, modifier=False):
    deltametric_vector_derivative = np.zeros((d + 1, len(prop)))
    try:
        metal_ind = mol.findMetal()[0]
        w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
        deltametric_vector_derivative = func(mol, w, metal_ind, d, oct=oct,
                                             catoms=catoms)
    except:
        print('Error, no metal found in mol object!')
        return False
    return (deltametric_vector_derivative)


def metal_only_deltametric(mol, prop, d, oct=True, catoms=None,
                           func=deltametric, modifier=False):
    deltametric_vector = np.zeros(d + 1)
    try:
        metal_ind = mol.findMetal()[0]
        w = construct_property_vector(mol, prop, oct=oct, modifier=modifier)
        deltametric_vector = func(mol, w, metal_ind, d, oct=oct,
                                  catoms=catoms)
    except:
        print('Error, no metal found in mol object!')
        return False
    return (deltametric_vector)


# Get the ligand_misc_descriptors
# custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
##                                    ax_con_int_list ,eq_con_int_list
# @param mol, mol3D class
# @param loud, bool, print out statements for debugging
# @return results_dictionary, vector, ax vs eq. charge (from OBMol) and denticity
def generate_all_ligand_misc(mol, loud, custom_ligand_dict=False):
    result_ax = list()
    result_eq = list()
    colnames = ['dent', 'charge']
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
            ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
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
    result_ax_charge = False
    result_eq_charge = False
    # loop over axial ligands
    if n_ax > 0:
        for i in range(0, n_ax):
            if mol.bo_dict:
                ax_ligand_list[i].mol.convert2OBMol2()
            else:
                ax_ligand_list[i].mol.convert2OBMol()
            if not (i == 0):
                result_ax_dent += ax_ligand_list[i].dent
                result_ax_charge += ax_ligand_list[i].mol.OBMol.GetTotalCharge()
            else:
                result_ax_dent = ax_ligand_list[i].dent
                result_ax_charge = ax_ligand_list[i].mol.OBMol.GetTotalCharge()
        # average axial results
        result_ax_dent = np.divide(result_ax_dent, n_ax)
        result_ax_charge = np.divide(result_ax_charge, n_ax)
    # loop over eq ligands
    if n_eq > 0:
        for i in range(0, n_eq):
            if mol.bo_dict:
                eq_ligand_list[i].mol.convert2OBMol2()
            else:
                eq_ligand_list[i].mol.convert2OBMol()
            if not (i == 0):
                result_eq_dent += eq_ligand_list[i].dent
                result_eq_charge += eq_ligand_list[i].mol.OBMol.GetTotalCharge()
            else:
                result_eq_dent = eq_ligand_list[i].dent
                result_eq_charge = eq_ligand_list[i].mol.OBMol.GetTotalCharge()
        # average eq results
        result_eq_dent = np.divide(result_eq_dent, n_eq)
        result_eq_charge = np.divide(result_eq_charge, n_eq)
        # save the results
    result_ax.append(result_ax_dent)
    result_ax.append(result_ax_charge)
    result_eq.append(result_eq_dent)
    result_eq.append(result_eq_charge)
    results_dictionary = {'colnames': colnames,
                          'result_ax': result_ax, 'result_eq': result_eq}
    return results_dictionary


## this function takes a
## symmetric (axial == axial,
## equatorial == equatorial)
## octahedral complex
## and returns autocorrelations for
## the axial an equatorial ligands
## custom_ligand_dict allows the user to skip the breakdown
## in cases where 3D geo is not correct/formed
## custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
##                                    ax_con_int_list ,eq_con_int_list
## with types: eq/ax_ligand_list list of mol3D
##             eq/ax_con_int_list list of list/tuple of int e.g,  [[1,2] [1,2]]
# @param mol, mol3D class, octhadral complex
# @param loud, bool, print lots of statements
# @param depth, int, number of steps 
# @return results_dictionary, dict, ligand autocorrelations
def generate_all_ligand_autocorrelations(mol, loud, depth=4, name=False, flag_name=False,
                                         custom_ligand_dict=False, NumB=False, Zeff=False):
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+=['effective_nuclear_charge']
        labels_strings += ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    result_ax_full = list()
    result_eq_full = list()
    result_ax_con = list()
    result_eq_con = list()
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
            ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
            mol, liglist, ligdents, ligcons, loud, name=False)
    else:
        ax_ligand_list = custom_ligand_dict["ax_ligand_list"]
        eq_ligand_list = custom_ligand_dict["eq_ligand_list"]
        ax_con_int_list = custom_ligand_dict["ax_con_int_list"]
        eq_con_int_list = custom_ligand_dict["eq_con_int_list"]
    ## count ligands
    n_ax = len(ax_ligand_list)
    n_eq = len(eq_ligand_list)
    colnames = []
    for ii, properties in enumerate(allowed_strings):
        ############### replaced find_ligand_autocorrelations_oct function here
        ## get full ligand AC
        ax_ligand_ac_full = []
        eq_ligand_ac_full = []
        for i in range(0, n_ax):
            if not list(ax_ligand_ac_full):
                ax_ligand_ac_full = full_autocorrelation(ax_ligand_list[i].mol, properties, depth)
            else:
                ax_ligand_ac_full += full_autocorrelation(ax_ligand_list[i].mol, properties, depth)
        ax_ligand_ac_full = np.divide(ax_ligand_ac_full, n_ax)
        for i in range(0, n_eq):
            if not list(eq_ligand_ac_full):
                eq_ligand_ac_full = full_autocorrelation(eq_ligand_list[i].mol, properties, depth)
            else:
                eq_ligand_ac_full += full_autocorrelation(eq_ligand_list[i].mol, properties, depth)
        eq_ligand_ac_full = np.divide(eq_ligand_ac_full, n_eq)
        ax_ligand_ac_con = []
        eq_ligand_ac_con = []
        for i in range(0, n_ax):
            if not list(ax_ligand_ac_con):
                ax_ligand_ac_con = atom_only_autocorrelation(ax_ligand_list[i].mol, properties, depth, ax_con_int_list[i])
            else:
                ax_ligand_ac_con += atom_only_autocorrelation(ax_ligand_list[i].mol, properties, depth, ax_con_int_list[i])
        ax_ligand_ac_con = np.divide(ax_ligand_ac_con, n_ax)
        for i in range(0, n_eq):
            if not list(eq_ligand_ac_con):
                eq_ligand_ac_con = atom_only_autocorrelation(eq_ligand_list[i].mol, properties, depth, eq_con_int_list[i])
            else:
                eq_ligand_ac_con += atom_only_autocorrelation(eq_ligand_list[i].mol, properties, depth, eq_con_int_list[i])
        eq_ligand_ac_con = np.divide(eq_ligand_ac_con, n_eq)
        ################
        this_colnames = []
        for i in range(0, depth + 1):
            this_colnames.append(labels_strings[ii] + '-' + str(i))
        colnames.append(this_colnames)
        result_ax_full.append(ax_ligand_ac_full)
        result_eq_full.append(eq_ligand_ac_full)
        result_ax_con.append(ax_ligand_ac_con)
        result_eq_con.append(eq_ligand_ac_con)
    if flag_name:
        results_dictionary = {'colnames': colnames, 'result_ax_full_ac': result_ax_full,
                              'result_eq_full_ac': result_eq_full,
                              'result_ax_con_ac': result_ax_con, 'result_eq_con_ac': result_eq_con}
    else:
        results_dictionary = {'colnames': colnames, 'result_ax_full': result_ax_full, 'result_eq_full': result_eq_full,
                              'result_ax_con': result_ax_con, 'result_eq_con': result_eq_con}
    return results_dictionary


## this function takes a
## symmetric (axial == axial,
## equatorial == equatorial)
## octahedral complex
## and returns autocorrelation derivatives for
## the axial an equatorial ligands
## custom_ligand_dict allows the user to skip the breakdown
## in cases where 3D geo is not correct/formed
## custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
##                                    ax_con_int_list ,eq_con_int_list
## with types: eq/ax_ligand_list list of mol3D
##             eq/ax_con_int_list list of list/tuple of int e.g,  [[1,2] [1,2]]
# @param mol, mol3D class, octhadral complex
# @param loud, bool, print lots of statements
# @param depth, int, number of steps 
# @return results_dictionary, dict, ligand autocorrelations
def generate_all_ligand_autocorrelation_derivatives(mol, loud, depth=4, name=False, flag_name=False,
                                                    custom_ligand_dict=False, NumB=False, Zeff=False):
    result_ax_full = None
    result_eq_full = None
    result_ax_con = None
    result_eq_con = None
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings += ['effective_nuclear_charge']
        labels_strings += ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
            ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
            mol, liglist, ligdents, ligcons, loud, name=False)
    else:
        ax_ligand_list = custom_ligand_dict["ax_ligand_list"]
        eq_ligand_list = custom_ligand_dict["eq_ligand_list"]
        ax_con_int_list = custom_ligand_dict["ax_con_int_list"]
        eq_con_int_list = custom_ligand_dict["eq_con_int_list"]
    ## count ligands
    n_ax = len(ax_ligand_list)
    n_eq = len(eq_ligand_list)
    for ii, properties in enumerate(allowed_strings):
        ## allocate the full jacobian matrix
        ax_full_j = np.zeros([depth + 1, mol.natoms])
        eq_full_j = np.zeros([depth + 1, mol.natoms])
        ax_con_j = np.zeros([depth + 1, mol.natoms])
        eq_con_j = np.zeros([depth + 1, mol.natoms])
        #################
        # full ligand ACs
        for i in range(0, n_ax):  # for each ax ligand
            ax_ligand_ac_full_derivative = full_autocorrelation_derivative(ax_ligand_list[i].mol, properties, depth)
            ## now we need to map back to full positions
            for jj, row in enumerate(ax_ligand_ac_full_derivative):
                for original_ids in list(ax_ligand_list[i].ext_int_dict.keys()):
                    ax_full_j[jj, original_ids] += np.divide(row[ax_ligand_list[i].ext_int_dict[original_ids]], n_ax)
        for i in range(0, n_eq):  # for each eq ligand
            ## now we need to map back to full positions
            eq_ligand_eq_full_derivative = full_autocorrelation_derivative(eq_ligand_list[i].mol, properties, depth)
            for jj, row in enumerate(eq_ligand_eq_full_derivative):
                for original_ids in list(eq_ligand_list[i].ext_int_dict.keys()):
                    eq_full_j[jj, original_ids] += np.divide(row[eq_ligand_list[i].ext_int_dict[original_ids]], n_eq)
        # ligand connection ACs
        for i in range(0, n_ax):
            ax_ligand_ac_con_derivative = atom_only_autocorrelation_derivative(ax_ligand_list[i].mol, properties, depth,
                                                                            ax_con_int_list[i])
            ## now we need to map back to full positions
            for jj, row in enumerate(ax_ligand_ac_con_derivative):
                for original_ids in list(ax_ligand_list[i].ext_int_dict.keys()):
                    ax_con_j[jj, original_ids] += np.divide(row[ax_ligand_list[i].ext_int_dict[original_ids]], n_ax)
        for i in range(0, n_eq):
            eq_ligand_ac_con_derivative = atom_only_autocorrelation_derivative(eq_ligand_list[i].mol, properties, depth,
                                                                            eq_con_int_list[i])
            ## now we need to map back to full positions
            for jj, row in enumerate(eq_ligand_ac_con_derivative):
                for original_ids in list(eq_ligand_list[i].ext_int_dict.keys()):
                    eq_con_j[jj, original_ids] += np.divide(row[eq_ligand_list[i].ext_int_dict[original_ids]], n_eq)
        ax_ligand_ac_full, eq_ligand_ac_full, ax_ligand_ac_con, eq_ligand_ac_con = ax_full_j, eq_full_j, ax_con_j, eq_con_j
        #################
        for i in range(0, depth + 1):
            colnames.append(['d' + labels_strings[ii] + '-' + str(i) + '/d' + labels_strings[ii] + str(j) for j in
                             range(0, mol.natoms)])
        if result_ax_full is None:
            result_ax_full = ax_ligand_ac_full
        else:
            result_ax_full = np.row_stack([result_ax_full, ax_ligand_ac_full])

        if result_eq_full is None:
            result_eq_full = eq_ligand_ac_full
        else:
            result_eq_full = np.row_stack([result_eq_full, eq_ligand_ac_full])

        if result_ax_con is None:
            result_ax_con = ax_ligand_ac_con
        else:
            result_ax_con = np.row_stack([result_ax_con, ax_ligand_ac_con])

        if result_eq_con is None:
            result_eq_con = eq_ligand_ac_con
        else:
            result_eq_con = np.row_stack([result_eq_con, eq_ligand_ac_con])
    if flag_name:
        results_dictionary = {'colnames': colnames, 'result_ax_full_ac': result_ax_full,
                              'result_eq_full_ac': result_eq_full,
                              'result_ax_con_ac': result_ax_con, 'result_eq_con_ac': result_eq_con}
    else:
        results_dictionary = {'colnames': colnames, 'result_ax_full': result_ax_full, 'result_eq_full': result_eq_full,
                              'result_ax_con': result_ax_con, 'result_eq_con': result_eq_con}
    return results_dictionary

## this function takes a
## symmetric (axial == axial,
## equatorial == equatorial)
## octahedral complex
## and returns deltametric autocorrelations for
## the axial an equatorial ligands
## custom_ligand_dict allows the user to skip the breakdown
## in cases where 3D geo is not correct/formed
## custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
##                                    ax_con_int_list ,eq_con_int_list
## with types: eq/ax_ligand_list list of mol3D
##             eq/ax_con_int_list list of list/tuple of int e.g,  [[1,2] [1,2]]
# @param mol, mol3D class, octhadral complex
# @param loud, bool, print lots of statements
# @param depth, int, number of steps 
# @return results_dictionary, dict, ligand autocorrelations
def generate_all_ligand_deltametrics(mol, loud, depth=4, name=False, flag_name=False,
                                     custom_ligand_dict=False, NumB=False, Zeff=False):
    result_ax_con = list()
    result_eq_con = list()
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
            ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
            mol, liglist, ligdents, ligcons, loud, name=False)
    else:
        ax_ligand_list = custom_ligand_dict["ax_ligand_list"]
        eq_ligand_list = custom_ligand_dict["eq_ligand_list"]
        ax_con_int_list = custom_ligand_dict["ax_con_int_list"]
        eq_con_int_list = custom_ligand_dict["eq_con_int_list"]
    ## count ligands
    n_ax = len(ax_ligand_list)
    n_eq = len(eq_ligand_list)
    for ii, properties in enumerate(allowed_strings):
        ####################
        ## get partial ligand AC
        ax_ligand_ac_con = []
        eq_ligand_ac_con = []
        for i in range(0, n_ax):
            if not list(ax_ligand_ac_con):
                ax_ligand_ac_con = atom_only_deltametric(ax_ligand_list[i].mol, properties, depth, ax_con_int_list[i])
            else:
                ax_ligand_ac_con += atom_only_deltametric(ax_ligand_list[i].mol, properties, depth, ax_con_int_list[i])
        ax_ligand_ac_con = np.divide(ax_ligand_ac_con, n_ax)
        for i in range(0, n_eq):
            if not list(eq_ligand_ac_con):
                eq_ligand_ac_con = atom_only_deltametric(eq_ligand_list[i].mol, properties, depth, eq_con_int_list[i])
            else:
                eq_ligand_ac_con += atom_only_deltametric(eq_ligand_list[i].mol, properties, depth, eq_con_int_list[i])
        eq_ligand_ac_con = np.divide(eq_ligand_ac_con, n_eq)
        ####################
        this_colnames = []
        for i in range(0, depth + 1):
            this_colnames.append(labels_strings[ii] + '-' + str(i))
        colnames.append(this_colnames)
        result_ax_con.append(ax_ligand_ac_con)
        result_eq_con.append(eq_ligand_ac_con)
    if flag_name:
        results_dictionary = {'colnames': colnames, 'result_ax_con_del': result_ax_con,
                              'result_eq_con_del': result_eq_con}
    else:
        results_dictionary = {'colnames': colnames, 'result_ax_con': result_ax_con, 'result_eq_con': result_eq_con}
    return results_dictionary


def generate_all_ligand_deltametric_derivatives(mol, loud, depth=4, name=False, flag_name=False,
                                                custom_ligand_dict=False, NumB=False, Zeff=False):
    ## custom_ligand_dict.keys() must be eq_ligands_list, ax_ligand_list
    ##                                    ax_con_int_list ,eq_con_int_list
    ## with types: eq/ax_ligand_list list of mol3D
    ##             eq/ax_con_int_list list of list/tuple of int e.g,  [[1,2] [1,2]]

    result_ax_con = None
    result_eq_con = None
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    if not custom_ligand_dict:
        liglist, ligdents, ligcons = ligand_breakdown(mol)
        ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, \
            ax_con_list, eq_con_list, built_ligand_list = ligand_assign(
            mol, liglist, ligdents, ligcons, loud, name=False)
    else:
        ax_ligand_list = custom_ligand_dict["ax_ligand_list"]
        eq_ligand_list = custom_ligand_dict["eq_ligand_list"]
        ax_con_int_list = custom_ligand_dict["ax_con_int_list"]
        eq_con_int_list = custom_ligand_dict["eq_con_int_list"]
    ## count ligands
    n_ax = len(ax_ligand_list)
    n_eq = len(eq_ligand_list)
    for ii, properties in enumerate(allowed_strings):
        ## allocate the full jacobian matrix
        ax_con_j = np.zeros([depth + 1, mol.natoms])
        eq_con_j = np.zeros([depth + 1, mol.natoms])
        #################
        for i in range(0, n_ax):
            ax_ligand_ac_con_derivative = atom_only_deltametric_derivative(ax_ligand_list[i].mol, properties, depth,
                                                                        ax_con_int_list[i])
            for jj, row in enumerate(ax_ligand_ac_con_derivative):
                for original_ids in list(ax_ligand_list[i].ext_int_dict.keys()):
                    ax_con_j[jj, original_ids] += np.divide(row[ax_ligand_list[i].ext_int_dict[original_ids]], n_ax)
        for i in range(0, n_eq):
            eq_ligand_ac_con_derivative = atom_only_deltametric_derivative(eq_ligand_list[i].mol, properties, depth,
                                                                        eq_con_int_list[i])
            for jj, row in enumerate(eq_ligand_ac_con_derivative):
                for original_ids in list(eq_ligand_list[i].ext_int_dict.keys()):
                    eq_con_j[jj, original_ids] += np.divide(row[eq_ligand_list[i].ext_int_dict[original_ids]], n_eq)
        #################
        ax_ligand_ac_con, eq_ligand_ac_con = ax_con_j, eq_con_j
        for i in range(0, depth + 1):
            colnames.append(['d' + labels_strings[ii] + '-' + str(i) + '/d' + labels_strings[ii] + str(j) for j in
                             range(0, mol.natoms)])
        if result_ax_con is None:
            result_ax_con = ax_ligand_ac_con
        else:
            result_ax_con = np.row_stack([result_ax_con, ax_ligand_ac_con])
        if result_eq_con is None:
            result_eq_con = eq_ligand_ac_con
        else:
            result_eq_con = np.row_stack([result_eq_con, eq_ligand_ac_con])
    if flag_name:
        results_dictionary = {'colnames': colnames, 'result_ax_con_del': result_ax_con,
                              'result_eq_con_del': result_eq_con}
    else:
        results_dictionary = {'colnames': colnames, 'result_ax_con': result_ax_con, 'result_eq_con': result_eq_con}
    return results_dictionary


def generate_metal_autocorrelations(mol, loud, depth=4, oct=True, flag_name=False,
                                    modifier=False, NumB=False, Zeff=False):
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = list()
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    for ii, properties in enumerate(allowed_strings):
        metal_ac = metal_only_autocorrelation(mol, properties, depth, oct=oct, modifier=modifier)
        this_colnames = []
        for i in range(0, depth + 1):
            this_colnames.append(labels_strings[ii] + '-' + str(i))
        colnames.append(this_colnames)
        result.append(metal_ac)
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_mc_ac': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_autocorrelation_derivatives(mol, loud, depth=4, oct=True, flag_name=False,
                                               modifier=False, NumB=False, Zeff=False):
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = None
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    for ii, properties in enumerate(allowed_strings):
        metal_ac_der = metal_only_autocorrelation_derivative(mol, properties, depth, oct=oct, modifier=modifier)
        for i in range(0, depth + 1):
            colnames.append(['d' + labels_strings[ii] + '-' + str(i) + '/d' + labels_strings[ii] + str(j) for j in
                             range(0, mol.natoms)])

        if result is None:
            result = metal_ac_der
        else:
            result = np.row_stack([result, metal_ac_der])
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_mc_ac': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_deltametrics(mol, loud, depth=4, oct=True, flag_name=False,
                                modifier=False, NumB=False, Zeff=False):
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = list()
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings+= ['effective_nuclear_charge']
        labels_strings+= ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    for ii, properties in enumerate(allowed_strings):
        metal_ac = metal_only_deltametric(mol, properties, depth, oct=oct, modifier=modifier)
        this_colnames = []
        for i in range(0, depth + 1):
            this_colnames.append(labels_strings[ii] + '-' + str(i))
        colnames.append(this_colnames)
        result.append(metal_ac)
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_mc_del': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_deltametric_derivatives(mol, loud, depth=4, oct=True, flag_name=False,
                                           modifier=False, NumB=False, Zeff=False):
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = None
    colnames = []
    allowed_strings = ['electronegativity', 'nuclear_charge', 'ident', 'topology', 'size']
    labels_strings = ['chi', 'Z', 'I', 'T', 'S']
    if Zeff:
        allowed_strings += ['effective_nuclear_charge']
        labels_strings += ['Zeff']
    if NumB:
        allowed_strings += ["num_bonds"]
        labels_strings += ["NumB"]
    for ii, properties in enumerate(allowed_strings):
        metal_ac_der = metal_only_deltametric_derivative(mol, properties, depth, oct=oct, modifier=modifier)
        for i in range(0, depth + 1):
            colnames.append(['d' + labels_strings[ii] + '-' + str(i) + '/d' + labels_strings[ii] + str(j) for j in
                             range(0, mol.natoms)])
        if result is None:
            result = metal_ac_der
        else:
            result = np.row_stack([result, metal_ac_der])
    if flag_name:
        results_dictionary = {'colnames': colnames, 'results_mc_del': result}
    else:
        results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary

################## Needed?

def generate_metal_ox_autocorrelations(oxmodifier, mol, loud, depth=4, oct=True, flag_name=False):
    ## oxmodifier - dict, used to modify prop vector (e.g. for adding
    ##             ONLY used with  ox_nuclear_charge    ox or charge)
    ##              {"Fe":2, "Co": 3} etc, normally only 1 metal...
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = list()
    colnames = []
    metal_ox_ac = metal_only_autocorrelation(mol, 'ox_nuclear_charge', depth, oct=oct, modifier=oxmodifier)
    this_colnames = []
    for i in range(0, depth + 1):
        this_colnames.append('O' + '-' + str(i))
    colnames.append(this_colnames)
    result.append(metal_ox_ac)
    results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_ox_autocorrelation_derivatives(oxmodifier, mol, loud, depth=4, oct=True, flag_name=False):
    ## oxmodifier - dict, used to modify prop vector (e.g. for adding 
    ##             ONLY used with  ox_nuclear_charge    ox or charge)
    ##              {"Fe":2, "Co": 3} etc, normally only 1 metal... 
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = None
    colnames = []
    metal_ox_ac = metal_only_autocorrelation_derivative(mol, 'ox_nuclear_charge', depth, oct=oct, modifier=oxmodifier)
    for i in range(0, depth + 1):
        colnames.append(['d' + 'O' + '-' + str(i) + '/d' + 'O' + str(j) for j in range(0, mol.natoms)])
    result = metal_ox_ac
    results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_ox_deltametrics(oxmodifier, mol, loud, depth=4, oct=True, flag_name=False):
    ## oxmodifier - dict, used to modify prop vector (e.g. for adding
    ##             ONLY used with  ox_nuclear_charge    ox or charge)
    ##              {"Fe":2, "Co": 3} etc, normally only 1 metal...
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = list()
    colnames = []
    metal_ox_ac = metal_only_deltametric(mol, 'ox_nuclear_charge', depth, oct=oct, modifier=oxmodifier)
    this_colnames = []
    for i in range(0, depth + 1):
        this_colnames.append('O' + '-' + str(i))
    colnames.append(this_colnames)
    result.append(metal_ox_ac)
    results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary


def generate_metal_ox_deltametric_derivatives(oxmodifier, mol, loud, depth=4, oct=True, flag_name=False):
    ## oxmodifier - dict, used to modify prop vector (e.g. for adding 
    ##             ONLY used with  ox_nuclear_charge    ox or charge)
    ##              {"Fe":2, "Co": 3} etc, normally only 1 metal... 
    #	oct - bool, if complex is octahedral, will use better bond checks
    result = list()
    colnames = []
    metal_ox_ac = metal_only_deltametric_derivative(mol, 'ox_nuclear_charge', depth, oct=oct, modifier=oxmodifier)
    for i in range(0, depth + 1):
        colnames.append(['d' + 'O' + '-' + str(i) + '/d' + 'O' + str(j) for j in range(0, mol.natoms)])

    result = metal_ox_ac
    results_dictionary = {'colnames': colnames, 'results': result}
    return results_dictionary