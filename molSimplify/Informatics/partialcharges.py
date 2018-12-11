from molSimplify.Classes.mol3D import *
import numpy as np
import glob

def fpriority(mol):
    # setting up variables
    fpriority_list = []
    fidx_list = []
    sidx_list = []
    satno_list = []
    ref_list = []
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

    for satno in sorted(set(satno_list)):
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
        if tidx_list == test_list:
            exit_signal = True

    for i in range(len(satno_list)):
        atno_list = []
        atno_list.append(str(satno_list[i]))
        if tatno_list[i] == []:
            atno_list.append('')
        else:
            atno_list.append(tatno_list[i][0])
        a = '.'.join(atno_list)
        fpriority_list.append(float(a))   
        
    return fpriority_list

def fsym(mol):
    # getting idxs of interest
    midx = mol.findMetal()[0] # monometallic complexes
    fidx_list = mol.getBondedAtoms(midx) # list of idx of the first-coord sphere
    fsym_list = []
    for idx in fidx_list:
        sym = mol.getAtom(idx).sym
        fsym_list.append(sym)
    
    return fsym_list    

def fvalency(mol):
    # getting idxs of interest
    midx = mol.findMetal()[0] # monometallic complexes
    fidx_list = mol.getBondedAtoms(midx) # list of idx of the first-coord sphere
    fvalency_list = []
    for idx in fidx_list:
        valency = len(mol.getBondedAtoms(idx)) - 1
        fvalency_list.append(valency)
    
    return fvalency_list    

def fcharge(mol,charge,bond=False):
    # getting idxs of interest
    midx = mol.findMetal()[0] # monometallic complexes
    fidx_list = mol.getBondedAtoms(midx) # list of idx of the first-coord sphere
    mol.calcCharges(charge,bond)
    fcharge_list = []
    for idx in fidx_list:
        partialq = mol.partialcharges[idx]
        fcharge_list.append(float(partialq))
    
    return fcharge_list

def scharge_ave(mol,charge,bond=False):
    mol.calcCharges(charge)
    # getting idxs of interest
    midx = mol.findMetal()[0] # monometallic complexes
    fidx_list = mol.getBondedAtoms(midx) # list of idx of the first-coord sphere
    mol.calcCharges(charge,bond)
    sidx_list = [mol.getBondedAtoms(fidx) for fidx in fidx_list]
    scharge_ave_list = []
    for i in range(len(sidx_list)):
        partialq = 0
        for j in range(len(sidx_list[i])):
            idx = sidx_list[i][j]
            if idx is not midx:
                partialq =+ mol.partialcharges[idx]
        charge_ave = partialq/len(sidx_list[i])
        scharge_ave_list.append(float(charge_ave))

    return scharge_ave_list

def fdistance(mol):
    # getting idxs of interest
    midx = mol.findMetal()[0] # monometallic complexes
    mcoord = mol.getAtom(midx).coords()
    fidx_list = mol.getBondedAtoms(midx) # list of idx of the first-coord sphere
    fdistance_list = []
    for idx in fidx_list:
        fcoord = mol.getAtom(idx).coords()
        d = distance(mcoord,fcoord)
        fdistance_list.append(float(d))
    
    return fdistance_list

def all_prop(mol,charge,bond=False):
    fprio_list = fpriority(mol)
    # fsym_list = fsym(mol)
    fva_list = fvalency(mol)
    fq_list = fcharge(mol,charge,bond)
    sq_ave_list = scharge_ave(mol,charge,bond)
    fd_list = fdistance(mol)
    prop_list = [fprio_list,fq_list,sq_ave_list,fva_list,fd_list]

    return prop_list

def f_prop(mol,charge,bond=False):
    fprio_list = fpriority(mol)
    # fsym_list = fsym(mol)
    fva_list = fvalency(mol)
    fq_list = fcharge(mol,charge,bond)
    sq_ave_list = scharge_ave(mol,charge,bond)
    fd_list = fdistance(mol)
    prop_list = [fprio_list,fq_list,fva_list]

    return prop_list

def features(mol,charge,bond=False):
    feature = []
    feature_names = []
    prop_list = all_prop(mol,charge,bond)
    midx = mol.findMetal()[0]
    manto = mol.getAtom(midx).atno
    a = np.array(prop_list)
    b = a.T[a[0].argsort()].T
    feature_list = b.tolist()
    feature_list[0] = [int(str(i).split('.')[0]) for i in feature_list[0]]
    feature.append(manto)
    for i in range(len(feature_list)):
        for j in range(len(feature_list[i])):
            feature.append(feature_list[i][j])
    feature_names = ['mato','fatno_1','fatno_2','fatno_3','fatno_4','fatno_5','fatno_6','fq_1','fq_2','fq_3','fq_4',
    'fq_5','fq_6','sqave_1','sqave_2','sqave_3','sqave_4','sqave_5','sqave_6','fval_1','fval_2','fval_3','fval_4',
    'fval_5','fval_6','bl_1','bl_2','bl_3','bl_4','bl_5','bl_6']

    return feature_names, feature
    # return feature

def ffeatures(mol,charge,bond=False):
    feature = []
    feature_names = []
    prop_list = f_prop(mol,charge,bond)
    midx = mol.findMetal()[0]
    manto = mol.getAtom(midx).atno
    a = np.array(prop_list)
    b = a.T[a[0].argsort()].T
    feature_list = b.tolist()
    feature_list[0] = [int(str(i).split('.')[0]) for i in feature_list[0]]
    feature.append(manto)
    for i in range(len(feature_list)):
        for j in range(len(feature_list[i])):
            feature.append(feature_list[i][j])
    feature_names = ['mato','fatno_1','fatno_2','fatno_3','fatno_4','fatno_5','fatno_6','fq_1','fq_2','fq_3','fq_4',
    'fq_5','fq_6','fval_1','fval_2','fval_3','fval_4','fval_5','fval_6']

    return feature_names, feature

