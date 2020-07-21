import numpy as np

from molSimplify.Classes.mol3D import mol3D
from molSimplify.Scripts.geometry import distance

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
            for sidx in xyz.getBondedAtoms(fidx):
                sidx_list.append([sidx])

    for i in range(len(fidx_list)):
        for fidx in fidx_list[i]:
            for j in range(len(sidx_list)):
                for sidx in sidx_list[j]:
                    BO = int(BOMatrix[fidx][sidx])
                    if BO == 0:
                        BO = 1
                    satno_str = str(xyz.getAtom(sidx).atno)
                    satno_list.append(int(BO * satno_str))

    for satno in sorted(set(satno_list)):
        satnocount = satno_list.count(satno)
        if satnocount > 1:
            s_sel_list = [i for i, atno in enumerate(
                satno_list) if atno is satno]
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
                for tidx in xyz.getBondedAtoms(sidx):
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
                sorted(ls, reverse=True)
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
    midx = xyz.findMetal()[0]  # monometallic complexes
    # list of idx of the first-coord sphere
    fidx_list = xyz.getBondedAtoms(midx)
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
    midx = xyz.findMetal()[0]  # monometallic complexes
    # list of idx of the first-coord sphere
    fidx_list = xyz.getBondedAtoms(midx)
    fvalency_list = []
    for idx in fidx_list:
        valency = len(xyz.getBondedAtoms(idx)) - 1
        fvalency_list.append(valency)

    return fvalency_list


def fcharge(xyzf, charge):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    xyz.calccharges(charge)
    # getting idxs of interest
    midx = xyz.findMetal()[0]  # monometallic complexes
    # list of idx of the first-coord sphere
    fidx_list = xyz.getBondedAtoms(midx)
    fcharge_list = []
    for idx in fidx_list:
        charge = xyz.partialcharges[idx]
        fcharge_list.append(float(charge))

    return fcharge_list


def scharge_ave(xyzf, charge):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    xyz.calccharges(charge)
    # getting idxs of interest
    midx = xyz.findMetal()[0]  # monometallic complexes
    # list of idx of the first-coord sphere
    fidx_list = xyz.getBondedAtoms(midx)
    sidx_list = [xyz.getBondedAtoms(fidx) for fidx in fidx_list]
    scharge_ave_list = []
    for i in range(len(sidx_list)):
        charge = 0
        for j in range(len(sidx_list[i])):
            idx = sidx_list[i][j]
            if idx is not midx:
                charge = + xyz.partialcharges[idx]
        charge_ave = charge/len(sidx_list[i])
        scharge_ave_list.append(float(charge_ave))

    return scharge_ave_list


def fdistance(xyzf):
    # setting properties
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    # getting idxs of interest
    midx = xyz.findMetal()[0]  # monometallic complexes
    mcoord = xyz.getAtom(midx).coords()
    # list of idx of the first-coord sphere
    fidx_list = xyz.getBondedAtoms(midx)
    fdistance_list = []
    for idx in fidx_list:
        fcoord = xyz.getAtom(idx).coords()
        d = distance(mcoord, fcoord)
        fdistance_list.append(float(d))

    return fdistance_list


def all_prop(xyzf, charge):
    fprio_list = fpriority(xyzf)
    # fsym_list = fsym(xyzf)
    fva_list = fvalency(xyzf)
    fq_list = fcharge(xyzf, charge)
    sq_ave_list = scharge_ave(xyzf, charge)
    fd_list = fdistance(xyzf)
    prop_list = [fprio_list, fq_list, sq_ave_list, fva_list, fd_list]

    return prop_list


def features(xyzf, charge):
    prop_list = all_prop(xyzf, charge)
    xyz = mol3D()
    xyz.readfromxyz(xyzf)
    midx = xyz.findMetal()[0]
    manto = xyz.getAtom(midx).atno
    a = np.array(prop_list)
    b = a.T[a[0].argsort()].T
    feature_list = b.tolist()
    feature_list[0] = [int(str(i).split('.')[0]) for i in feature_list[0]]
    feature = []
    feature.append(manto)
    for i in range(len(feature_list)):
        for j in range(len(feature_list[i])):
            feature.append(feature_list[i][j])

    return feature
