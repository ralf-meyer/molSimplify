from molSimplify.Classes.mol3D import *
import numpy as np
import glob

# xyzf = '/Users/tzuhsiungyang/Downloads/cuacetate1k2acetate1o-pyridylphenyl1_0_1_RIJCOSX-B3LYP-D3_BS_TS-ar-carboxylation.numfreq.xyz'

def fpriority(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    # setting up variables
    fpriority_list = []
    fidx_list = []
    sidx_list = []
    satno_list = []
    ref_list = []
    exit_signal = True
    # getting bond-order matrix
    xyz.convert2OBMol()
    BOMatrix = xyz.populateBOMatrix()

    # preping for the loop
    fidx_list.append(xyz.findMetal())
    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            for sidx in xyz.getBondedAtomsBOMatrix(fidx):
                sidx_list.append([sidx])

    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            for j in range(len(sidx_list)):
                for sidx in sidx_list[j]:
                    BO = int(BOMatrix[fidx][sidx])
                    satno_str = str(xyz.getAtom(sidx).atno)
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
                for tidx in xyz.getBondedAtomsBOMatrix(sidx):
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
                        tatno_str = str(xyz.getAtom(tidx).atno)
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

def fsym(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    # getting idxs of interest
    midx = xyz.findMetal()[0] # monometallic complexes
    fidx_list = xyz.getBondedAtomsBOMatrix(midx) # list of idx of the first-coord sphere
    fsym_list = []
    for idx in fidx_list:
        sym = xyz.getAtom(idx).sym
        fsym_list.append(sym)
    
    return fsym_list    

def fvalency(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    # getting idxs of interest
    midx = xyz.findMetal()[0] # monometallic complexes
    fidx_list = xyz.getBondedAtomsBOMatrix(midx) # list of idx of the first-coord sphere
    fvalency_list = []
    for idx in fidx_list:
        valency = len(xyz.getBondedAtomsBOMatrix(idx)) - 1
        fvalency_list.append(valency)
    
    return fvalency_list    

def fcharge(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    xyz.calccharges()
    # getting idxs of interest
    midx = xyz.findMetal()[0] # monometallic complexes
    fidx_list = xyz.getBondedAtomsBOMatrix(midx) # list of idx of the first-coord sphere
    fcharge_list = []
    for idx in fidx_list:
        charge = xyz.partialcharges[idx]
        fcharge_list.append(float(charge))
    
    return fcharge_list

def scharge_ave(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    xyz.calccharges()
    # getting idxs of interest
    midx = xyz.findMetal()[0] # monometallic complexes
    fidx_list = xyz.getBondedAtomsBOMatrix(midx) # list of idx of the first-coord sphere
    sidx_list = [xyz.getBondedAtomsBOMatrix(fidx) for fidx in fidx_list]
    scharge_ave_list = []
    for i in range(len(sidx_list)):
        charge = 0
        for j in range(len(sidx_list[i])):
            idx = sidx_list[i][j]
            if idx is not midx:
                charge =+ xyz.partialcharges[idx]
        charge_ave = charge/len(sidx_list[i])
        scharge_ave_list.append(float(charge_ave))

    return scharge_ave_list

def fdistance(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    # getting idxs of interest
    midx = xyz.findMetal()[0] # monometallic complexes
    mcoord = xyz.getAtom(midx).coords()
    fidx_list = xyz.getBondedAtomsBOMatrix(midx) # list of idx of the first-coord sphere
    fdistance_list = []
    for idx in fidx_list:
        fcoord = xyz.getAtom(idx).coords()
        d = distance(mcoord,fcoord)
        fdistance_list.append(float(d))
    
    return fdistance_list

def all_prop(xyzf):
    fprio_list = fpriority(xyzf)
    fsym_list = fsym(xyzf)
    fva_list = fvalency(xyzf)
    fq_list = fcharge(xyzf)
    sq_ave_list = scharge_ave(xyzf)
    fd_list = fdistance(xyzf)
    prop_list = [fprio_list,fq_list,sq_ave_list,fva_list,fd_list]

    return prop_list

def features(xyzf):
    prop_list = all_prop(xyzf)
    a = np.array([])
    i_size = len(prop_list)
    j_size = len(prop_list[0])
    for i in range(i_size):
        b = np.array(prop_list[i])
        a = np.append(a,b)
        b = np.reshape(a,(i_size,j_size))
        c = b.T[a[:j_size].argsort()].T
        feature_list = c.tolist()
        feature_list[0] = [int(str(i).split('.')[0]) for i in feature_list[0]]

    return feature_list
