## @file nn_prep.py
#  Helper routines for ANN integration
#  
#   Written by Nick Ynag for HJK Group
#  
#  Dpt of Chemical Engineering, MIT

from molSimplify.Classes import mol3D
from molSimplify.Informatics.autocorrelation import *
from molSimplify.Informatics.graph_analyze import *
from molSimplify.Informatics.partialcharges import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.globalvars import *
from sklearn.kernel_ridge import KernelRidge
import numpy as np
import pandas as pd
import time, re
from sets import Set

## wrapper to get KRR predictions for bondl_core3D, bondl_m3D, bondl_m3Dsub from a known mol3D using partial charges
#  @param mol mol3D of the molecule
#  @param charge charge of the molecule
#  @return KRR-predicted bondl_core3D
#  KRR accuracies for bondl_core3D: 98.2% (training score) and 47.6 (test score)
#  KRR accuracies for bondl_m3D: 99.5% (training score) and 51.1 (test score)
def invoke_KRR_from_mol3d_dQ(mol,charge):
    # # find the metal from RACs 
    # metal = mol.getAtom(mol.findMetal()[0]).symbol()
    # ox_modifier = {metal:oxidation_state}
    # get partialQs
    feature_names,features = ffeatures(mol,charge)
    # # get one-hot-encoding (OHE)
    # descriptor_names,descriptors = create_OHE(descriptor_names,descriptors, metal,oxidation_state)
    # # set exchange fraction
    # descriptor_names += ['alpha']
    # descriptors += [alpha]
    ## KRR initiation
    # defined variables
    globs = globalvars()
    if globs.custom_path: # test if a custom path is used:
         X_train_csv = str(globs.custom_path).rstrip('/') + "/python_krr/X_train_TS.csv"
         y_train_csv = str(globs.custom_path).rstrip('/') + "/python_krr/y_train_TS.csv"
    else:
        X_train_csv = resource_filename(Requirement.parse("molSimplify"),"molSimplify/python_krr/X_train_TS.csv")
        y_train_csv = resource_filename(Requirement.parse("molSimplify"),"molSimplify/python_krr/y_train_TS.csv")
    X_train = pd.read_csv(X_train_csv,header=None)
    y_train = pd.read_csv(y_train_csv,header=None)
    kernel = 'rbf'
    keys = []
    bondls = []
    for targets in ['bondl_core3D','bondl_m3D']:#,'bondl_m3Dsub']:
        keys.append(targets)
        if targets == 'bondl_core3D':
        # KRR parameters for bondl_core3D 
            alpha = 0.1
            gamma = 4.6415888336127775
            mean_y_train = 1.8556069976566096
            std_y_train = 0.08511267085380758
            mean_X_train = np.array([1.1886128903870394, 1.0746595698697274, 1.0089390403652372, 1.0051636435711488, 
                                    0.9639844597149281, 1.5924309727104378])
            std_X_train = np.array([1.4887238067607071, 1.4391120341824508, 1.351343230273359, 1.302911028297482,
                                    1.1511093513567663, 0.7366350688359029])

        if targets == 'bondl_m3D':
        # KRR parameters for bondl_core3D 
            alpha = 0.015848931924611134
            gamma = 8.531678524172808
            mean_y_train = 1.1429284052746633
            std_y_train = 0.04763054722349127
            mean_X_train = np.array([-1.17136495, -1.09058534, -1.04062806, -1.01379334, -0.92612448, -1.30558513])
            std_X_train = np.array([1.36359461, 1.32785945, 1.26392399, 1.21494676, 1.0253893, 0.5940198])

        # model initation
        X_test = np.array(features[7:13])
        X_test = (X_test - mean_X_train) / std_X_train
        model = KernelRidge(kernel=kernel,alpha=alpha,gamma=gamma)
        model.fit(X_train,y_train)
        y_test = model.predict([X_test])
        y_test = y_test * std_y_train + mean_y_train
        bondl = y_test[0][0]
        bondls.append(bondl)
    
    bondl_dict = dict(zip(keys,bondls))

    return bondl_dict

## wrapper to get KRR predictions for bondl_core3D from a known mol3D using RAC-190
#  @param mol mol3D of the molecule
#  @param charge charge of the molecule
#  @return KRR-predicted bondl_core3D
#  KRR accuracies: 98.2% (training score) and 47.6 (test score)
def invoke_KRR_from_mol3d_RACs(mol,charge):
    # # find the metal from RACs 
    # metal = mol.getAtom(mol.findMetal()[0]).symbol()
    # ox_modifier = {metal:oxidation_state}
    # get partialQs
    feature_names,features = ffeatures(mol,charge)
    # # get one-hot-encoding (OHE)
    # descriptor_names,descriptors = create_OHE(descriptor_names,descriptors, metal,oxidation_state)
    # # set exchange fraction
    # descriptor_names += ['alpha']
    # descriptors += [alpha]
    ## KRR initiation
    # defined variables
    X_train = pd.read_csv('/Users/tzuhsiungyang/anaconda2/envs/molSimplify/molSimplify/molSimplify/python_krr/X_train_TS.csv',header=None)
    y_train = pd.read_csv('/Users/tzuhsiungyang/anaconda2/envs/molSimplify/molSimplify/molSimplify/python_krr/y_train_TS.csv',header=None)
    kernel = 'rbf'
    alpha = 0.1
    gamma = 4.6415888336127775
    mean_y_train = 1.8556069976566096
    std_y_train = 0.08511267085380758
    mean_X_train = np.array([1.1886128903870394, 1.0746595698697274, 1.0089390403652372, 1.0051636435711488, 0.9639844597149281, 1.5924309727104378])
    std_X_train = np.array([1.4887238067607071, 1.4391120341824508, 1.351343230273359, 1.302911028297482, 1.1511093513567663, 0.7366350688359029])
    # model initation
    X_test = np.array(features[7:13])
    X_test = (X_test - mean_X_train) / std_X_train
    model = KernelRidge(kernel=kernel,alpha=alpha,gamma=gamma)
    model.fit(X_train,y_train)
    y_test = model.predict([X_test])
    y_test = y_test * std_y_train + mean_y_train
    bondl_core3D = y_test[0][0]

    return bondl_core3D

## Gets the RACs of a given atidx
#  @param mol mol3D of this molecule
#  @param atidx the index of the atom of concern
#  @return descriptor_names updated names
#  @return descriptors updated RACs
def get_descriptor_vector_for_atidx(mol, atidx, depth=4):
    descriptor_names = []
    descriptors = []
    result_dictionary = generate_atomonly_autocorrelations(mol, atidx,False, depth)
    for colnames in result_dictionary['colnames']:
        descriptor_names += colnames
    for results in result_dictionary['results']:
        descriptors += results.tolist()
    result_dictionary = generate_atomonly_deltametrics(mol, atidx, False, depth)
    for colnames in result_dictionary['colnames']:
        descriptor_names += colnames
    for results in result_dictionary['results']:
        descriptors += results.tolist()

    return descriptor_names, descriptors

