## @file nn_prep.py
#  Helper routines for ANN integration
#  
#   Written by Nick Yang for HJK Group
#  
#  Dpt of Chemical Engineering, MIT

from molSimplify.Classes import mol3D
from molSimplify.Informatics.autocorrelation import *
from molSimplify.Informatics.graph_analyze import *
from molSimplify.Informatics.partialcharges import *
from molSimplify.Classes.mol3D import *
from molSimplify.Classes.globalvars import *
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.kernel_ridge import KernelRidge
from sklearn.multioutput import MultiOutputRegressor
import numpy as np
import csv, glob, os
import matplotlib.pyplot as plt

def feature_prep(mol, idx):
    # setting up variables
    fidx_list = []
    sidx_list = []
    satno_list = []
    ref_list = []
    fd_list = []
    fa_list = []
    idx_list = [0] * 6
    exit_signal = True
    # getting bond-order matrix
    mol.convert2OBMol()
    BOMatrix = mol.populateBOMatrix()

    # preping for the loop
    fidx_list.append(mol.findMetal())
    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            for sidx in mol.getBondedAtoms(fidx):
                sidx_list.append([sidx])

    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            for j in range(len(sidx_list)):
                for sidx in sidx_list[j]:
                    BO = int(BOMatrix[fidx][sidx])
                    if BO == 0:
                        BO = 1
                    satno_str = str(mol.getAtom(sidx).atno)
                    satno_list.append(int(BO * satno_str))

    for satno in set(satno_list):
        satnocount = satno_list.count(satno)
        if satnocount > 1:
            s_sel_list = [i for i,atno in enumerate(satno_list) if atno is satno]
            exit_signal = False

    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            ref_list.append(fidx)

    # starting the loop
    tidx_list = []
    tatno_list = []
    for i in range(len(sidx_list)):
        tidx_list.append([])
        tatno_list.append([])

    while not exit_signal:
        fpriority_list = []
        for i in s_sel_list:
            t_list = []
            for sidx in sidx_list[i]:
                for tidx in mol.getBondedAtoms(sidx):
                    if tidx not in ref_list:
                        t_list.append(tidx)
            tidx_list[i] = t_list
        # print(sidx_list)
        # print(tidx_list)
        for i in s_sel_list:
            for sidx in sidx_list[i]:
                atno_list = tatno_list[i]
                ls = []
                for j in s_sel_list:
                    for tidx in tidx_list[j]:
                        BO = int(BOMatrix[sidx][tidx])
                        tatno_str = str(mol.getAtom(tidx).atno)
                        ls.append(BO * tatno_str)
                sorted(ls,reverse=True)
                for j in ls:
                    atno_list.append(j)
                a = ''.join(atno_list)
            tatno_list[i] = [a]
        sidx_list = []
        for i in range(len(tidx_list)):
            sidx_list.append(tidx_list[i])
        for i in s_sel_list:
            for sidx in sidx_list[i]:
                ref_list.append(sidx)
        test_list = []
        for i in range(len(sidx_list)):
            test_list.append([])
        # get priorities
        for i in range(len(satno_list)):
            atno_list = []
            atno_list.append(str(satno_list[i]))
            if tatno_list[i] == []:
                atno_list.append('')
            else:
                atno_list.append(tatno_list[i][0])
            a = '.'.join(atno_list)
            fpriority_list.append(float(a))
        if tidx_list == test_list or len(set(fpriority_list)) == 6:
        # if tidx_list == test_list:
            exit_signal = True
    # get distance
    # idx = np.argsort(np.array(fpriority_list))[-1]
    sidx_list = mol.getBondedAtomsByCoordNo(fidx_list[0][0], 6)
    refcoord = mol.getAtom(sidx_list[idx]).coords()
    mcoord = mol.getAtom(fidx_list[0][0]).coords()
    idx0 = 0
    dist5 = 0
    idx5 = 0
    idx1_4 = []
    fprio1_4 = []
    sxyzs = []
    ssd_list = []
    for i, sidx in enumerate(sidx_list):
        sxyz = mol.getAtom(sidx).coords()
        dist = distance(refcoord,sxyz)
        if i == idx:
            idx0 = i
        elif dist > dist5:
            dist5 = dist
            idx5 = i
        idx1_4.append(i)
        fprio1_4.append(fpriority_list[i])
        sxyzs.append(sxyz)
        ssd_list.append(dist)
        fd_list.append(distance(mcoord,sxyz))
    idx1_4.pop(idx0)
    idx1_4.pop(idx5)
    fprio1_4.pop(idx0)
    fprio1_4.pop(idx5)
    idx_list[0] = idx0
    idx_list[5] = idx5
    idx1 = idx1_4[np.argsort(np.array(fprio1_4))[3]]
    sxyz1 = sxyzs[idx1]
    idx2_ = idx1_4[np.argsort(np.array(fprio1_4))[2]]
    sxyz2_ = sxyzs[idx2_]
    idx3_ = idx1_4[np.argsort(np.array(fprio1_4))[1]]
    sxyz3_ = sxyzs[idx3_]
    idx4_ = idx1_4[np.argsort(np.array(fprio1_4))[0]]
    sxyz4_ = sxyzs[idx4_]
    fd1_4 = []
    fd1_4.append(distance(sxyz1, sxyz1))
    fd1_4.append(distance(sxyz1, sxyz2_))
    fd1_4.append(distance(sxyz1, sxyz3_))
    fd1_4.append(distance(sxyz1, sxyz4_))
    idx3 = idx1_4[np.argsort(np.array(fd1_4))[-1]] + idx1 - 4
    if idx3 == idx2_:
        if fpriority_list[idx3_] > fpriority_list[idx4_]:
            idx2 = idx3_
            idx4 = idx4_
        else:
            idx2 = idx4_
            idx4 = idx2_
    elif idx3 == idx4_:
        if fpriority_list[idx2_] > fpriority_list[idx3_]:
            idx2 = idx2_
            idx4 = idx3_
        else:
            idx2 = idx3_
            idx4 = idx2_
    else:
        if fpriority_list[idx2_] > fpriority_list[idx4_]:
            idx2 = idx2_
            idx4 = idx4_
        else:
            idx2 = idx4_
            idx4 = idx2_
    # get ax, eq, ax idxes
    idx_list[1] = idx1
    idx_list[2] = idx2
    idx_list[3] = idx3
    idx_list[4] = idx4
    fpriority_list = np.array(fpriority_list)[idx_list].tolist()
    fd_list = np.array(fd_list)[idx_list].tolist()

    return fpriority_list, fd_list, idx_list

def normalize(data, mean, std):
    data_norm = (data - mean) / std
    data_norm = np.nan_to_num(data_norm)

    return data_norm

## predict labels using krr with a given csv file
#  @param csvf the csv file containing headers (first row), data, and label
#  @param colnum_label the column number for the label column
#  @param colnum_desc the starting column number for the descriptor columns
#  @return y_train_data, y_train_pred, y_test_data, y_test_pred, score
def ML_prediction_krr(csvf, colnum_label, colnum_desc):
    # read in desc and label
    f = open(csvf, 'r')
    fcsv = csv.reader(f)
    headers = np.array(next(f, None).rstrip('\r\n').split(','))
    X = []
    y = []
    for line in fcsv:
        if len(line) == 309:
            descs = []
            for desc in line[colnum_desc:]:
                descs.append(float(desc))
            X.append(descs)
            y.append(float(line[colnum_label]))
    X = np.array(X)
    y = np.array(y)
    ## process desc and label
    mean_X = np.mean(X, axis=0)
    std_X = np.std(X, axis=0)
    mean_y = np.mean(y, axis=0)
    std_y = np.std(y, axis=0)
    X_norm = normalize(X, mean_X, std_X)
    y_norm = normalize(y, mean_y, std_y)
    # split to train and test
    X_norm_train, X_norm_test, y_norm_train, y_norm_test = train_test_split(X_norm, y_norm, test_size=0.2, random_state=0)
    ## end
    # feature selection
    selector = RandomForestRegressor(random_state=0, n_estimators=100)
    selector.fit(X_norm_train, y_norm_train)
    X_norm_train_impts = selector.feature_importances_
    idxes = np.where(X_norm_train_impts > 0.01)[0]
    features_sel = headers[idxes]
    X_norm_train_sel = X_norm_train.T[idxes].T
    X_norm_test_sel = X_norm_test.T[idxes].T
    ## training with krr
    # krr parameters
    kernel = 'rbf'
    gamma = 1
    alpha = 1
    factor_lower = 0.5
    factor_higher = 2
    gamma_lower = gamma * factor_lower
    gamma_higher = gamma * factor_higher
    alpha_lower = alpha * factor_lower
    alpha_higher = alpha * factor_higher
    lin = 9
    # optimize hyperparameters
    while gamma == 1 or alpha == 1 or \
            (gamma < gammas[lin / 2 - 1] or gamma > gammas[lin / 2]) or \
            (alpha < alphas[lin / 2 - 1] or alpha > alphas[lin / 2]):
        gammas = np.linspace(gamma_lower, gamma_higher, lin)
        alphas = np.linspace(alpha_lower, alpha_higher, lin)
        tuned_parameters = [{'kernel': [kernel], 'gamma': gammas, 'alpha': alphas}]
        regr = GridSearchCV(KernelRidge(), tuned_parameters, cv=5, scoring='neg_mean_absolute_error')
        regr.fit(X_norm_train_sel, y_norm_train)
        gamma = regr.best_params_['gamma']
        alpha = regr.best_params_['alpha']
        # factor_lower *= 2
        # factor_higher *= 0.5
        gamma_lower = gamma * factor_lower
        gamma_higher = gamma * factor_higher
        alpha_lower = alpha * factor_lower
        alpha_higher = alpha * factor_higher
    # final model
    regr = KernelRidge(kernel=kernel, alpha=alpha, gamma=gamma)
    regr.fit(X_norm_train_sel, y_norm_train)
    # predictions
    y_norm_train_pred = regr.predict(X_norm_train_sel)
    y_train_pred = y_norm_train_pred * std_y + mean_y
    y_train_data = y_norm_train * std_y + mean_y
    y_norm_test_pred = regr.predict(X_norm_test_sel)
    y_test_pred = y_norm_test_pred * std_y + mean_y
    y_test_data = y_norm_test * std_y + mean_y
    # performance
    score_train = regr.score(X_norm_train_sel, y_norm_train)
    score_test = regr.score(X_norm_test_sel, y_norm_test)
    MAE_train = mean_absolute_error(y_train_data, y_train_pred)
    MAE_test = mean_absolute_error(y_test_data, y_test_pred)
    stat_names = ['score_train', 'score_test', 'MAE_train', 'MAE_test']
    stats = [score_train, score_test, MAE_train, MAE_test]
    stat_dict = dict(zip(stat_names, stats))

    return y_train_data, y_train_pred, y_test_data, y_test_pred, stat_dict

## predict labels using gradient boosting regressor (GBR) with a given csv file
#  @param csvf the csv file containing headers (first row), data, and label
#  @param colnum_i_label the starting column number for the label column
#  @param colnum_j_label the ending column number for the label column + 1
#  @param colnum_desc the starting column number for the descriptor columns
#  @return y_train_data, y_train_pred, y_test_data, y_test_pred, score
def ML_prediction_gbr(csvf, colnum_i_label, colnum_j_label, colnum_desc):
    # read in desc and label
    f = open(csvf, 'r')
    fcsv = csv.reader(f)
    headers = np.array(next(f, None).rstrip('\r\n').split(','))
    X = []
    y = []
    for line in fcsv:
        if len(line) == 309:
            descs = []
            labels = []
            for desc in line[colnum_desc:]:
                descs.append(float(desc))
            for label in line[colnum_i_label:colnum_j_label]:
                labels.append(float(label))
            X.append(descs)
            y.append(labels)
    X = np.array(X)
    y = np.array(y)
    ## process desc and label
    mean_X = np.mean(X, axis=0)
    std_X = np.std(X, axis=0)
    mean_y = np.mean(y, axis=0)
    std_y = np.std(y, axis=0)
    X_norm = normalize(X, mean_X, std_X)
    y_norm = normalize(y, mean_y, std_y)
    # split to train and test
    X_norm_train, X_norm_test, y_norm_train, y_norm_test = train_test_split(X_norm, y_norm, test_size=0.2, random_state=0)
    ## end
    # feature selection
    # selector = RandomForestRegressor(random_state=0, n_estimators=100)
    # selector.fit(X_norm_train, y_norm_train)
    # X_norm_train_impts = selector.feature_importances_
    # idxes = np.where(X_norm_train_impts > 0.01)[0]
    # features_sel = headers[idxes]
    idxes = range(len(X_norm_train.T))
    X_norm_train_sel = X_norm_train.T[idxes].T
    X_norm_test_sel = X_norm_test.T[idxes].T
    ## training with gbr
    regr = MultiOutputRegressor(GradientBoostingRegressor(random_state=0))
    # final model
    regr.fit(X_norm_train_sel, y_norm_train)
    # predictions
    y_norm_train_pred = regr.predict(X_norm_train_sel)
    y_train_pred = y_norm_train_pred * std_y + mean_y
    y_train_data = y_norm_train * std_y + mean_y
    y_norm_test_pred = regr.predict(X_norm_test_sel)
    y_test_pred = y_norm_test_pred * std_y + mean_y
    y_test_data = y_norm_test * std_y + mean_y
    # performance
    score_train = regr.score(X_norm_train_sel, y_norm_train)
    score_test = regr.score(X_norm_test_sel, y_norm_test)
    MAE_train = mean_absolute_error(y_train_data, y_train_pred)
    MAE_test = mean_absolute_error(y_test_data, y_test_pred)
    stat_names = ['score_train', 'score_test', 'MAE_train', 'MAE_test']
    stats = [score_train, score_test, MAE_train, MAE_test]
    stat_dict = dict(zip(stat_names, stats))

    return y_train_data, y_train_pred, y_test_data, y_test_pred, stat_dict

## wrapper to get KRR predictions for bondl_core3D, bondl_m3D, bondl_m3Dsub from a known mol3D using partial charges
#  @param mol mol3D of the molecule
#  @param charge charge of the molecule
#  @return KRR-predicted bondl_core3D
#  KRR accuracies for bondl_core3D: 98.2% (training score) and 47.6 (test score)
#  KRR accuracies for bondl_m3D: 99.5% (training score) and 51.1 (test score)
def invoke_KRR_from_mol3d_dQ(mol,charge):
    X_norm_train = []
    y_norm_train = []
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
         X_norm_train_csv = str(globs.custom_path).rstrip('/') + "/python_krr/X_norm_train_TS.csv"
         y_norm_train_csv = str(globs.custom_path).rstrip('/') + "/python_krr/y_norm_train_TS.csv"
    else:
        X_norm_train_csv = resource_filename(Requirement.parse("molSimplify"),"molSimplify/python_krr/X_norm_train_TS.csv")
        y_norm_train_csv = resource_filename(Requirement.parse("molSimplify"),"molSimplify/python_krr/y_norm_train_TS.csv")
    f = open(X_norm_train_csv,'r')
    for line in csv.reader(f):
        X_norm_train.append([float(i) for i in line])
    X_norm_train = np.array(X_norm_train)
    f = open(y_norm_train_csv,'r')
    for line in csv.reader(f):
        y_norm_train.append([float(i) for i in line])
    y_norm_train = np.array(y_norm_train)
    # X_norm_train = pd.read_csv(X_norm_train_csv,header=None)
    # y_norm_train = pd.read_csv(y_norm_train_csv,header=None)
    kernel = 'rbf'
    keys = []
    bondls = []
    for targets in ['bondl_core3D','bondl_m3D']:#,'bondl_m3Dsub']:
        keys.append(targets)
        if targets == 'bondl_core3D':
        # KRR parameters for bondl_core3D 
            alpha = 0.1
            gamma = 4.6415888336127775
            mean_y_norm_train = 1.8556069976566096
            std_y_norm_train = 0.08511267085380758
            mean_X_norm_train = np.array([1.1886128903870394, 1.0746595698697274, 1.0089390403652372, 1.0051636435711488, 
                                    0.9639844597149281, 1.5924309727104378])
            std_X_norm_train = np.array([1.4887238067607071, 1.4391120341824508, 1.351343230273359, 1.302911028297482,
                                    1.1511093513567663, 0.7366350688359029])

        if targets == 'bondl_m3D':
        # KRR parameters for bondl_core3D 
            alpha = 0.015848931924611134
            gamma = 8.531678524172808
            mean_y_norm_train = 1.1429284052746633
            std_y_norm_train = 0.04763054722349127
            mean_X_norm_train = np.array([-1.17136495, -1.09058534, -1.04062806, -1.01379334, -0.92612448, -1.30558513])
            std_X_norm_train = np.array([1.36359461, 1.32785945, 1.26392399, 1.21494676, 1.0253893, 0.5940198])

        # model initation
        X_norm_test = np.array(features[7:13])
        X_norm_test = (X_norm_test - mean_X_norm_train) / std_X_norm_train
        model = KernelRidge(kernel=kernel,alpha=alpha,gamma=gamma)
        model.fit(X_norm_train,y_norm_train)
        y_norm_test = model.predict([X_norm_test])
        y_norm_test = y_norm_test * std_y_norm_train + mean_y_norm_train
        bondl = y_norm_test[0][0]
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
    X_norm_train = pd.read_csv('/Users/tzuhsiungyang/anaconda2/envs/molSimplify/molSimplify/molSimplify/python_krr/X_norm_train_TS.csv',header=None)
    y_norm_train = pd.read_csv('/Users/tzuhsiungyang/anaconda2/envs/molSimplify/molSimplify/molSimplify/python_krr/y_norm_train_TS.csv',header=None)
    kernel = 'rbf'
    alpha = 0.1
    gamma = 4.6415888336127775
    mean_y_norm_train = 1.8556069976566096
    std_y_norm_train = 0.08511267085380758
    mean_X_norm_train = np.array([1.1886128903870394, 1.0746595698697274, 1.0089390403652372, 1.0051636435711488, 0.9639844597149281, 1.5924309727104378])
    std_X_norm_train = np.array([1.4887238067607071, 1.4391120341824508, 1.351343230273359, 1.302911028297482, 1.1511093513567663, 0.7366350688359029])
    # model initation
    X_norm_test = np.array(features[7:13])
    X_norm_test = (X_norm_test - mean_X_norm_train) / std_X_norm_train
    model = KernelRidge(kernel=kernel,alpha=alpha,gamma=gamma)
    model.fit(X_norm_train,y_norm_train)
    y_norm_test = model.predict([X_norm_test])
    y_norm_test = y_norm_test * std_y_norm_train + mean_y_norm_train
    bondl_core3D = y_norm_test[0][0]

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

def default_plot(x, y):
    # defs for plt
    # xlabel = r'% HF exchange'
    # ylabel = r'${\rm \Delta}{\rm E}^{\rm X-HS}$'
    colors = ['r', 'g', 'b', '.75', 'orange', 'k']
    markers = ['o', 's', 'D', 'v', '^', '<', '>']
    font = {'family': 'sans-serif',
            # 'weight' : 'bold',
            'size': 22}
    # figure size
    plt.figure(figsize=(7.5, 6))
    # dealing with axes
    # x = sorted(x)
    # y = sorted(y)
    # x_min = round(x[0],2)
    # x_max = round(x[-1],2)
    # plt.xlim(x_min, x_max)
    # y_min = round(y[0],2)
    # y_max = round(y[-1],2)
    # plt.ylim(y_min, y_max)
    # plt.xlabel(xlabel)
    # plt.ylabel(ylabel)
    # dealing with ticks
    ax = plt.axes()
    # ax.xaxis.set_major_locator(ticker.MultipleLocator((x_max - x_min) / 4))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator((x_max - x_min) / 8))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator((y_max - y_min) / 4))
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator((y_max - y_min) / 8))
    plt.tick_params(which='both', axis='both', direction='in', bottom=True, top=True, right=True, left=True)
    plt.rcParams['axes.linewidth'] = 3
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.major.width'] = 3
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.major.width'] = 3
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['xtick.minor.width'] = 3
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['ytick.minor.width'] = 3
    plt.tight_layout()

    plt.rc('font', **font)
    plt.plot(x, y, 'o', markeredgecolor='k')
    # plt.plot([x_min, x_max], [x_min, x_max], 'k', linestyle='dashed')
    plt.show()
# plt.imshow(data,interpolation='none')
# # plt.imshow(data,interpolation='nearest')
# plt.savefig('relative_energies_for_Fe-py4.eps',dpi=400)



