# from molSimplify.Classes.mol3D import *
# sub1 = mol3D()
# sub2 = mol3D()
# subtot = mol3D()
# sub1.readfromxyz('phenyl_-1.xyz')
# sub2.readfromxyz('phenyl_-1.xyz')
# sub1.translate([2,2,2])
# subtot.copymol3D(sub1)
# subtot.combine(sub2,[(5,16,1)])
# subtot.convert2OBMol()
# OBMol = subtot.OBMol
# ff = openbabel.OBForceField.FindForceField('mmff94')
# ff.Setup(OBMol)
# ff.SteepestDescent(50000)
# ff.GetCoordinates(OBMol)
# subtot.OBMol = OBMol
# subtot.convert2mol3D()
# subtot.writexyz('phenyl-phenyl.xyz')


# from molSimplify.Classes.partialcharges import *

# xyzf = glob.glob('ACNCFE_modified.RES1.xyz')[0]
# # setting properties
# xyz = mol3D()
# xyz.readfromxyz(xyzf)
# # setting up variables
# fpriority_list = []
# fidx_list = []
# sidx_list = []
# satno_list = []
# ref_list = []
# exit_signal = True
# # getting bond-order matrix
# xyz.convert2OBMol()
# BOMatrix = xyz.populateBOMatrix()

# # preping for the loop
# fidx_list.append(xyz.findMetal())
# for i in range(len(fidx_list)):
#     for fidx in fidx_list[i]:
#         for sidx in xyz.getBondedAtoms(fidx):
#             sidx_list.append([sidx])

# for i in range(len(fidx_list)):
#     for fidx in fidx_list[i]:
#         for j in range(len(sidx_list)):
#             for sidx in sidx_list[j]:
#                 BO = int(BOMatrix[fidx][sidx])
#                 if BO == 0:
#                     BO = 1
#                 satno_str = str(xyz.getAtom(sidx).atno)
#                 satno_list.append(int(BO * satno_str))

# for satno in sorted(set(satno_list)):
#     satnocount = satno_list.count(satno)
#     if satnocount > 1:
#         s_sel_list = [i for i,atno in enumerate(satno_list) if atno is satno]
#         exit_signal = False

# for i in range(len(fidx_list)):
#     for fidx in fidx_list[i]:
#         ref_list.append(fidx)

# # starting the loop
# tidx_list = []
# tatno_list = []
# for i in range(len(sidx_list)):
#     tidx_list.append([])
#     tatno_list.append([])

# for i in s_sel_list:
#     t_list = []
#     for sidx in sidx_list[i]:
#         for tidx in xyz.getBondedAtoms(sidx):
#             if tidx not in ref_list:
#                 t_list.append(tidx)
#     tidx_list[i] = t_list

# # print(sidx_list)
# # print(tidx_list)
# for i in s_sel_list:
#     for sidx in sidx_list[i]:
#         atno_list = tatno_list[i]
#         ls = []
#         for j in s_sel_list:
#             for tidx in tidx_list[j]:
#                 BO = int(BOMatrix[sidx][tidx])
#                 tatno_str = str(xyz.getAtom(tidx).atno)
#                 ls.append(BO * tatno_str)
#         sorted(ls,reverse=True)
#         for j in ls:
#             atno_list.append(j)
#         a = ''.join(atno_list)
#     tatno_list[i] = [a]


# with open('chenru_features_include_Q2.txt','r') as fin:
#     for i in fin.readlines():
#         feach = open('temp','w')
#         feach.write(i)
#         feach.close()
#         fop = open('temp','rb')
#         a = pickle.load(fop)
#         print(a)

# from molSimplify.Classes.mol3D import *
# import glob

# xyzf = '/Users/tzuhsiungyang/Runs/cu_2_en_1_formamidate_1_bromobenzene_multiple_cu_0_s_2/cu_2_en_1_formamidate_1_bromobenzene_0_11_cu_0_s_2/cu_2_en_1_formamidate_1_bromobenzene_0_11_cu_0_s_2.xyz'

# xyz = mol3D()
# xyz.readfromxyz(xyzf)
# midx_list = xyz.findMetal()
# fidx_list = []
# subidx_list = []

# for midx in midx_list:
#     for fidx in xyz.getBondedAtoms(midx):
#         if xyz.getAtom(fidx).sym != 'C' and xyz.getAtom(fidx).sym != 'Br':
#             fidx_list.append(fidx)
#         else:
#             subidx_list.append(fidx)


# for subidx in subidx_list:
#     subcoord = xyz.getAtom(subidx).coords()
#     f = [0,0,0]
#     for fidx in fidx_list:
#         f[0] += xyz.getAtom(fidx).coords()[0]/3
#         f[1] += xyz.getAtom(fidx).coords()[1]/3
#         f[2] += xyz.getAtom(fidx).coords()[2]/3
#     angle = vecangle(subcoord,f)
#     print('The angle between ' + str(subidx) + ' and centroid is ' + str(angle))

# import ast
# import numpy as np

# f = open('chenru_features_include_Q2.txt','r')
# a = np.array([])
# j = 0
# for line in f.readlines():
#     b_ls = ast.literal_eval(line)
#     i = len(b_ls)
#     j += 1
#     # for i in range(len(b_ls)):
#     #     b = np.asarray(b_ls[i])
#     #     a = np.append(a,b)
#     b = np.asarray(b_ls)
#     a = np.concatenate((a,b))

# np.reshape(a,(j,7))

# print(a)

# import ast
# import numpy as np

# f = open('chenru_features_include_Q2.txt','r')
# a = np.array([])
# j = 0
# len_list = []
# for line in f.readlines():
#     b_ls = ast.literal_eval(line)
#     len_list.append(len(b_ls))

# sorted(len_list,reverse=True)

## extract hessian from .hess and then diagonalize

import glob, os, re, sys
import numpy as np

# partial hessian

# extract hessian from .hess
hess = sys.argv[1]
f = open(hess,'r')
dof = int(f.readlines()[13].rstrip('\n'))
repeats = dof / 5 + 1
end_line = repeats * (dof + 1)
hessian = []
f = open(hess,'r')
for line in f.readlines()[14:end_line+14]:
    new_line = line.rstrip('\n').split(' ')
    H = [float(i) for i in new_line if i != '']
    hessian.append(H)

ls = []
for i in range(dof):
    ls.append([])

npf = np.asarray(ls)
for i in range(repeats):
    a = []
    for j in hessian[1+i*(dof+1):(i+1)*(dof+1)]:
        a.append(j[1:6])
    npa = np.asarray(a)
    npf = np.concatenate((npf,npa),axis=1)

idx_list = [11,13,12]
dof = []
for idx in idx_list:
    for i in range(idx*3,idx*3+3):
        dof.append(i)

partial_hess = []

for i in range(len(dof)):
    idx_i = dof[i]
    partial_hess.append([npf[idx][ii] for ii in dof])

# partial_hess = [[0.041677998743, 0.041179646923, 0.011033922837, 0.0012780504069, 0.0023894512887, 0.0043170271304, -0.026985523574, -0.043959595491, -0.015853455845],
#     [0.041179646923, 0.1116341609, 0.0078639300402, 0.0059849224541, -0.015229505848, 0.004878834599, -0.052726610788, -0.074365135965, -0.029028252688], 
#     [0.011033922837, 0.0078639300402, 0.09599160515, 0.027114763354, 0.027599925377, 0.00099830785582, -0.041304880141, -0.054636509255, -0.017319539258], 
#     [0.0012780504069, 0.0059849224541, 0.027114763354, -0.0069679064636, -0.024128860561, -0.026611331339, 0.0152634557, 0.033455851278, 0.01940617514], 
#     [0.0023894512887, -0.015229505848, 0.027599925377, -0.024128860561, -0.010867956478, -0.033130689958, 0.045632795266, 0.051813268426, 0.032639257168], 
#     [0.0043170271304, 0.004878834599, 0.00099830785582, -0.026611331339, -0.033130689958, 0.018384216132, 0.037284187187, 0.0503594692, 0.0010729231416], 
#     [-0.026985523574, -0.052726610788, -0.041304880141, 0.0152634557, 0.045632795266, 0.037284187187, 0.43100707883, -0.19705412973, -0.072901745246], 
#     [-0.043959595491, -0.074365135965, -0.054636509255, 0.033455851278, 0.051813268426, 0.0503594692, -0.19705412973, 0.27257520281, -0.095167291695], 
#     [-0.015853455845, -0.029028252688, -0.017319539258, 0.01940617514, 0.032639257168, 0.0010729231416, -0.072901745246, -0.095167291695, 0.51996030059]]

# extract hessian from .hess
hess = sys.argv[2]
f = open(hess,'r')
dof = int(f.readlines()[13].rstrip('\n'))
repeats = dof / 5 + 1
end_line = repeats * (dof + 1)
hessian = []
f = open(hess,'r')
for line in f.readlines()[14:end_line+14]:
    new_line = line.rstrip('\n').split(' ')
    H = [float(i) for i in new_line if i != '']
    hessian.append(H)

ls = []
for i in range(dof):
    ls.append([])

npf = np.asarray(ls)
for i in range(repeats):
    a = []
    for j in hessian[1+i*(dof+1):(i+1)*(dof+1)]:
        a.append(j[1:6])
    npa = np.asarray(a)
    npf = np.concatenate((npf,npa),axis=1)

# # get hessian elements for OHC
# idx_list = [6,7,8]
# dof = [range(idx*3,idx*3+3) for idx in idx_list]
# ls_all = []
# for i in range(len(dof)):
#     for idx_i in dof[i]:
#         ls = []
#         for j in range(len(dof)):
#             for idx_j in dof[j]:
#                 ls.append(npf[idx_i,idx_j])
#         ls_all.append(ls)

# replace selected elements of the hessian
# idx_list = [11,13,12]
# dof = [range(idx*3,idx*3+3) for idx in idx_list]
# for i in range(len(dof)):
#     for idx_i in dof[i]:
#         for j in range(len(dof)):
#             for idx_j in dof[j]:
#                 npf[idx_i,idx_j] = partial_hess[]
idx_list = [11,13,12]
dof = []
for idx in idx_list:
    for i in range(idx*3,idx*3+3):
        dof.append(i)

for i in range(len(dof)):
    idx_i = dof[i]
    for j in range(len(dof)):
        idx_j = dof[j]
        npf[idx_i,idx_j] = partial_hess[i][j]


# write the top portion of the old hessian file into the new one
f = open('try.hess','ab')
fold = open(hess,'r')
for i in fold.readlines()[:14]:
    f.write(i)

f = open('try.hess','ab')
# reshape the hessian for Orca
ncol = npf.shape[1]
for col in range(ncol/5):
    a = npf[col*5:col*5+5].T
    b = np.array([range(col*5,col*5+5)])
    a = np.concatenate((b,a),axis=0)
    b = np.array([[0]])
    c = np.array([range(ncol)]).T
    b = np.concatenate((b,c))
    a = np.concatenate((b,a),axis=1)
    a = a.astype(str)
    a[0,0] = ''
    a[0] = np.char.replace(a[0],'.0','')
    for i in range(ncol+1):
        a[i][0] = np.char.replace(a[i][0],'.0','')
    np.savetxt(f, a, delimiter=' ', fmt='%s')

f.close()

f = open('try.hess','ab')
if (ncol % 5) != 0:
    a = npf[ncol-(ncol % 5):ncol].T
    b = np.array([range(ncol-(ncol % 5),ncol)])
    a = np.concatenate((b,a),axis=0)
    b = np.array([[0]])
    c = np.array([range(ncol)]).T
    b = np.concatenate((b,c))
    a = np.concatenate((b,a),axis=1)
    a = a.astype(str)
    a[0,0] = ''
    a[0] = np.char.replace(a[0],'.0','')
    for i in range(ncol+1):
        a[i][0] = np.char.replace(a[i][0],'.0','')
    np.savetxt(f, a, delimiter=' ', fmt='%s')

f.close()

# write the bottom portion of the old hessian file into the new one
f = open('try.hess','r')
nlines_f = len(f.readlines())
f.close()
fold = open(hess,'r')
nlines_fold = len(fold.readlines())
fold.close()
f = open('try.hess','ab')
fold = open(hess,'r')
for i in fold.readlines()[nlines_f:nlines_fold]:
    f.write(i)

f.close()
fold.close()     

# ## read in cdxml
# from molSimplify.Scripts.io import *
# import glob, pybel
# cdxml = glob.glob('PMe3.cdxml')[0]
# obconv = openbabel.OBConversion() # ob Class
# obmol = openbabel.OBMol() # ob Class
# obconv.SetInFormat('cdxml') # ob Method to set cdxml
# obconv.ReadFile(obmol,cdxml) # ob Method to reaad cdxml into a OBMol()
# obmol.NumAtoms()
# idx_list = []
# atno_list = []
# for idx in range(obmol.NumAtoms()):
#     if obmol.GetAtom(idx+1).GetAtomicNum() == 26:
#         idx_list.append(idx)
#         atno_list.append(obmol.GetAtom(idx+1).GetAtomicNum())
#         obmol.GetAtom(idx+1).SetAtomicNum(14)

# pymol = pybel.Molecule(obmol)
# pymol.make3D()
# pymol.localopt()
# for i in range(len(idx_list)):
#     idx = idx_list[i]
#     atno = atno_list[i]
#     obmol.GetAtom(idx+1).SetAtomicNum(atno)

# fname = re.sub(r'.cdxml','.xyz',cdxml) # file name for the new xyz
# obconv.WriteFile(obmol,fname)

# # read in cdx and output xyz
# from molSimplify.Scripts.io import *
# import glob, pybel
# cdx = glob.glob('toluene.cdx')[0]
# obconv = openbabel.OBConversion() # ob Class
# obmol = openbabel.OBMol() # ob Class
# obconv.SetInFormat('cdx') # ob Method to set cdx
# obconv.ReadFile(obmol,cdx) # ob Method to reaad cdx into a OBMol()
# obmol.NumAtoms()
# idx_list = []
# atno_list = []
# for idx in range(obmol.NumAtoms()):
#     if obmol.GetAtom(idx+1).GetAtomicNum() == 26:
#         idx_list.append(idx)
#         atno_list.append(obmol.GetAtom(idx+1).GetAtomicNum())
#         obmol.GetAtom(idx+1).SetAtomicNum(14)

# pymol = pybel.Molecule(obmol)
# pymol.make3D()
# pymol.localopt()
# for i in range(len(idx_list)):
#     idx = idx_list[i]
#     atno = atno_list[i]
#     obmol.GetAtom(idx+1).SetAtomicNum(atno)

# fname = re.sub(r'.cdx','.xyz',cdx) # file name for the new xyz
# obconv.WriteFile(obmol,fname)


# obmol.AddHydrogens()
# obff = openbabel.OBForceField.FindForceField('mmff94')
# obff.Setup(obmol)
# obff.SteepestDescent(1000)
# obff.GetCoordinates(OBMol)
# fname = re.sub(r'.cdxml','.xyz',cdxml) # file name for the new xyz
# obconv.WriteFile(obmol,fname)

# # read in cdx and output cdxml
# from molSimplify.Scripts.io import *
# import glob, pybel
# cdx = glob.glob('toluene.cdx')[0]
# obconv = openbabel.OBConversion() # ob Class
# obmol = openbabel.OBMol() # ob Class
# obconv.SetInFormats('cdx') # ob Method to set cdx
# obconv.ReadFile(obmol,cdx) # ob Method to reaad cdx into a OBMol()

# fname = re.sub(r'.cdx','.cdxml',cdx) # file name for the new xyz
# obconv.WriteFile(obmol,fname)

# from molSimplify.Classes.mol3D import *
# import glob

# xyzf = glob.glob('*triazole*')[0]
# xyz = mol3D()
# xyz.readfromxyz(xyzf)
# xyz.BCM_opt(0,4,2.35)

# ## CDM

# from molSimplify.Classes.mol3D import *
# import glob

# reacting_at_list = [21,26,32]

# xyzf = glob.glob('/Users/tzuhsiungyang/Runs/pd_2_N-quinolinylbutyramidate_1_e2acetate_1_s_3/pd_2_N-quinolinylbutyramidate_1_e2acetate_1_s_3/pd*xyz')[0]
# xyz = mol3D()
# xyz.readfromxyz(xyzf)
# midx_list = xyz.findMetal()
# constr_list = []
# for midx in midx_list:
#     for fidx in xyz.getBondedAtomsBOMatrix(midx):
#         if fidx not in reacting_at_list:
#             constr_list.append(fidx)

# xyz.convert2OBMol()
# OBMol = xyz.OBMol
# for bond in openbabel.OBMolBondIter(OBMol):
#     idx1 = bond.GetBeginAtomIdx()
#     idx2 = bond.GetEndAtomIdx()
#     if (idx1-1 == 32 and idx2-1 == 0) or (idx2-1 == 32 and idx1-1 == 0):
#         bond.SetBO(0)
#     if (idx1-1 == 21 and idx2-1 == 26) or (idx2-1 == 21 and idx1-1 == 26):
#         bond.SetBO(0)
#     if (idx1-1 == 31 and idx2-1 == 33) or (idx2-1 == 31 and idx1-1 == 33):
#         bond.SetBO(2)

# xyz.OBMol = OBMol
# xyz.convert2mol3D()
# xyz.convert2OBMol()
# OBMol = xyz.OBMol
# ff = openbabel.OBForceField.FindForceField('uff')
# constr = openbabel.OBFFConstraints()
# for fidx in constr_list:
#     constr.AddAtomConstraint(fidx+1)

# for midx in midx_list:
#     constr.AddAtomConstraint(midx+1)

# constr.AddDistanceConstraint(22,27,1.521) # bondl_sub
# constr.AddDistanceConstraint(27,33,1.284) # bondl_m3D
# constr.AddDistanceConstraint(1,33,2.839) # bondl_core3D
# constr.AddDistanceConstraint(1,22,2.202)
# constr.AddDistanceConstraint(1,33,27,35.1) # bangle_m3D
# s = ff.Setup(OBMol,constr)
# ff.SteepestDescent(5000)
# ff.GetCoordinates(OBMol)
# xyz.OBMol = OBMol
# xyz.convert2mol3D()
# xyz.writexyz('/Users/tzuhsiungyang/Runs/pd_2_N-quinolinylbutyramidate_1_e2acetate_1_s_3/pd_2_N-quinolinylbutyramidate_1_e2acetate_1_s_3/debug.xyz')

# ## 

# from molSimplify.Classes.mol3D import *
# import glob

# reacting_at_list = [31,42]

# xyzf = glob.glob('/Users/tzuhsiungyang/Runs/pd_2_pme3_2_phenyl_2_s_3/pd_2_pme3_2_phenyl_2_s_3/pd*xyz')[0]
# xyz = mol3D()
# xyz.readfromxyz(xyzf)
# midx_list = xyz.findMetal()
# constr_list = []
# for midx in midx_list:
#     constr_list.append(midx)
#     for fidx in xyz.getBondedAtomsBOMatrix(midx):
#         if fidx not in reacting_at_list:
#             constr_list.append(fidx)

# xyz.convert2OBMol()
# OBMol = xyz.OBMol
# for bond in openbabel.OBMolBondIter(OBMol):
#     idx1 = bond.GetBeginAtomIdx()
#     idx2 = bond.GetEndAtomIdx()
#     if (idx1-1 in constr_list) and (idx2-1 in constr_list):
#         bond.SetBO(0)

# xyz.OBMol = OBMol
# xyz.convert2mol3D()
# xyz.convert2OBMol()
# OBMol = xyz.OBMol
# ff = openbabel.OBForceField.FindForceField('uff')
# constr = openbabel.OBFFConstraints()
# for fidx in constr_list:
#     constr.AddAtomConstraint(fidx+1)

# constr.AddDistanceConstraint(32,43,1.9) # bondl_sub
# constr.AddDistanceConstraint(1,32,2.0) # bondl_m3D
# constr.AddDistanceConstraint(1,43,2.0) # bondl_core3D
# constr.AddDistanceConstraint(1,32,43,54.5) # bangle_m3D
# s = ff.Setup(OBMol,constr)
# ff.SteepestDescent(5000)
# ff.GetCoordinates(OBMol)
# xyz.OBMol = OBMol
# xyz.convert2mol3D()
# xyz.writexyz('/Users/tzuhsiungyang/Runs/pd_2_pme3_2_phenyl_2_s_3/pd_2_pme3_2_phenyl_2_s_3/debug.xyz')

# sidx_generator = lambda x : x if x in midx_list else core3D.getBondedAtoms(subcatoms_ext[0])[0]
# sidx = [sidx_generator(i) for i in subcores.keys()]

# from molSimplify.Classes.partialcharges import *

# xyzf = '/Users/tzuhsiungyang/Dropbox (MIT)/Kulik Group Project Sean/all_oct_Fe_components_062618/JUPRAN_modified.RES1.xyz'
# features(xyzf,0)
# method='QEq'

# # def calccharges(self,method='QEq'):
# self.convert2OBMol()
# self.
# charge = openbabel.OBChargeModel.FindType(method)
# charge.ComputeCharges(self.OBMol)
# self.partialcharges = charge.GetPartialCharges()

# # frag xyz generator

# from molSimplify.Classes.mol3D import *
# from molSimplify.Classes.ligand import *
# import glob, itertools, re

# for xyzf in glob.glob('*xyz'):
#     xyz = mol3D()
#     xyz.readfromxyz(xyzf)
#     midx = xyz.findMetal()[0]
#     liglist,ligdent,ligcons = ligand_breakdown(xyz)
#     idx_list = range(len(liglist))
#     frag_idx = 0
#     cfor_list = []
#     draws = len(liglist)-1
#     for subset in itertools.combinations(idx_list,draws):
#         unique = True
#         dent = 0
#         frag_list = []
#         for i in subset:
#             dent_i = ligdent[i]
#             dent += dent_i
#             frag_list.append(liglist[i])
#         if dent == 5:
#             frag_xyz = mol3D()
#             frag_xyz.addAtom(xyz.getAtom(midx))
#             for i in range(len(frag_list)):
#                 for j in range(len(frag_list[i])):
#                     aidx = frag_list[i][j]
#                     frag_xyz.addAtom(xyz.getAtom(aidx))
#             frag_xyz.convert2OBMol()
#             cfor = frag_xyz.OBMol.GetFormula()
#             if cfor not in cfor_list:
#                 frag_fname = re.sub(r'.xyz','_frag' + str(frag_idx),xyzf)
#                 frag_xyz.writexyz(frag_fname)
#                 frag_idx += 1
#                 cfor_list.append(cfor)
#     print(xyzf + ' is complete.')

# import openbabel, re

# xyzf = '/Users/tzuhsiungyang/Downloads/ACPICF_modified.RES1.xyz'
# obconv = openbabel.OBConversion()
# obmol = openbabel.OBMol()
# obconv.SetInFormat('xyz')
# obconv.ReadFile(obmol,xyzf)
# obmol.AddHydrogens()
# xyzf_new = re.sub(r'.xyz','_addedH.xyz',xyzf)
# obconv.WriteFile(obmol,xyzf_new)

# from ccdc import io
# from ccdc.molecule import Atom,Bond,Molecule
# import re, glob

# for mol2 in glob.glob('*mol2'):
#     mol = io.MoleculeReader(mol2)[0]
#     if 'Unknown' not in str(mol.bonds):
#         mol.add_hydrogens()
#         mol2_new = re.sub(r'.mol2','_addedH.mol2',mol2)
#         f = open(mol2_new,'w')
#         f.write(mol.to_string('mol2'))
#         f.close()
#         print(mol2 + ' is complete.')


# from ccdc import io
# import glob

# for cif in glob.glob('*cif'):
#     ref = cif.split('_')[0]
#     csd_reader = io.EntryReader('CSD')
#     entry = csd_reader.entry(ref)
#     pub = entry.publication
#     pub.doi
#     # journal = pub.journal
#     # j_full_name = journal.full_name
#     # j_full_name.publisher_name


# # multiple gaussian fit

# from scipy.optimize import curve_fit
# from scipy.signal import argrelextrema
# import numpy as np
# import matplotlib.pyplot as plt
# import re

# def func(x, *params):
#     y = np.zeros_like(x)
#     for i in range(0, len(params), 3):
#         ctr = params[i]
#         amp = params[i+1]
#         wid = params[i+2]
#         y = y + amp * np.exp(-((x - ctr)/wid)**2)
#     return y

# fname = 'Fe-N_BL'

# data = np.loadtxt('/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/' + fname + '.csv', delimiter=',', unpack=True)
# x, y = data

# maxima = argrelextrema(y,np.greater)[0].tolist()
# guess = []
# for i in range(2):
#     idx = maxima.pop()
#     guess.append(x[idx])
#     guess.append(y[idx])
#     guess.append((x[idx+1]-x[idx-1])/1.5)

# lower_bound = [i*0.5 for i in guess]
# upper_bound = [i*2 for i in guess]

# popt, pcov = curve_fit(func, x, y, p0=guess,bounds=(lower_bound,upper_bound))
# print popt

# eqn = 'y = a*exp(-1/2*((x-b)/c)^2)'
# first_g = 'a = ' + str(round(popt[0],2)) + '\n'
# first_g += 'b = ' + str(round(popt[1],0)) + '\n'
# first_g += 'c = ' + str(round(popt[2],3))
# second_g = 'a = ' + str(round(popt[3],2)) + '\n'
# second_g += 'b = ' + str(round(popt[4],0)) + '\n'
# second_g += 'c = ' + str(round(popt[5],3))

# # # Mn-N text box
# # plt.text(1.9,800,eqn,bbox=dict(alpha=0.5))
# # plt.text(2.0,380,second_g)
# # plt.text(2.15,600,first_g)

# # Fe-N textbox
# plt.text(1.65,2100,eqn,bbox=dict(alpha=0.5))
# plt.text(1.83,1500,second_g)
# plt.text(2.15,1100,first_g)

# # # Co-N textbox
# # plt.text(2.1,2800,eqn,bbox=dict(alpha=0.5))
# # plt.text(1.83,2200,second_g)
# # plt.text(2.12,1300,first_g)

# plt.plot(x, y)
# x_smooth = np.linspace(x[0],x[-1],num=100)
# fit = func(x_smooth, *popt)
# plt.plot(x_smooth, fit , 'r-')
# xlabel_ = fname.split('_')[0]
# plt.xlabel(xlabel_)
# plt.ylabel('Frequency')
# # x_lower_upper = popt[0]-popt[2]*3
# # x_fill = np.linspace(x[0],x_lower_upper,num=100)
# # y_fill = func(x_fill, *popt)
# # plt.fill_between(x_fill,y_fill,0,facecolor='r')
# # x_fill_list.pop()
# # x_upper_lower = popt[3]+popt[5]*3
# # x_fill2 = np.linspace(x_upper_lower,x[-1],num=100)
# # y_fill2 = func(x_fill2, *popt)
# # plt.fill_between(x_fill2,y_fill2,0,facecolor='r')

# plt.show()

# f = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_selected_BL.txt'
# ele_dic = {'Mn':25,'Fe':26,'Co':27}
# ele = fname.split('-')[0]
# matno = ele_dic[ele]
# fop = open(f,'r')
# fout = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_selected_BL_with_label.txt'
# for line in fop.readlines():
#     line_ = line.rstrip('\n')
#     line_list = line_.split(',')
#     ref = line_list[0]
#     matno_f = line_list[1]
#     bl = line_list[2]
#     ss = ''
#     if len(line_list) > 3:
#         ss = ',' + line_list[3]
#     if int(matno_f) == (matno):
#         if bl < x_lower_upper:
#             ss = ',ls'
#         if bl > x_upper_lower:
#             ss = ',hs'
#     foutop = open(fout,'a')
#     foutop.write(ref + ',' + matno_f + ',' + str(bl) + ss + '\n')
#     foutop.close()

# # Calculate charge of the fragments

# from molSimplify.Classes.mol3D import *
# import sys,glob

# xyzf = sys.argv[1]
# ref = xyzf.splot('.')[0] + '.' + xyzf.splot('.')[1]
# chargef = '/home/nickyang/oct_1TM/first_row_v3/addedH/xyz_natom_le_30/xyzf_charge.txt'
# charge_r = open(chargef,'r')
# for line in charge_r.readlines():
#     line_rstrip = line.rstrip('\n')
#     xyzff = line_rstrip.split(' ')[0]
#     if xyzff == xyzf:
#         charge = line_rstrip.split(' ')[1]

# xyz_r = open(xyzf,'r')
# xyzlist = xyz_r.readlines()
# for fragf in glob.glob(ref + '_frag*xyz'):
#     frag_r = open(frag,'r')
#     xyz3D = mol3D()
#     tmp = open('tmp.xyz','w')
#     natms = 0
#     for i,line in enumerate(frag_r.readlines()):
#         if i > 2 and line in xyzlist:
#             tmp.write(line)
#             natoms += 1
#     tmp.insert(0,'Title')
#     tmp.insert(0,natms)
#     xyz3D.readfromxyz('tmp.xyz')
#     xyz3D.

# # 
# import openbabel,glob

# fout = open('ele_count.txt','w')

# # xyzf = '/home/nickyang/oct_1TM/first_row_v3/addedH/xyz_natom_le_15/FULLAZ-RES2-addedH_0_5_RC-cat.xyz'
# for xyzf_ in glob.glob('*RC-cat.xyz'):
#     xyzf = xyzf_.rstrip('\n')
# # xyzf = xyzf_.rstrip('\n')
#     ftype = xyzf.split('.')[-1]
#     ls = xyzf.split('_')
#     for i in range(3):
#         charge = ls.pop()
    
#     obconv = openbabel.OBConversion()
#     obmol = openbabel.OBMol()
#     obconv.SetInFormat(ftype)
#     obconv.ReadFile(obmol,xyzf)
#     obmol.SetTotalCharge(int(charge))
#     ele = obmol.GetTotalSpinMultiplicity()
#     fout.write(str(ele) + ',' + xyzf + '\n')



# # formula writer
# import openbabel,glob

# fout = open('xyzf_formula.txt','w')

# for xyzf_ in glob.glob('*xyz'):
#     # xyzf = '/home/nickyang/oct_1TM/first_row_v3/addedH/ADOSIW_modified.RES1_addedH.xyz'
#     xyzf = xyzf_.rstrip('\n')
#     ftype = xyzf.split('.')[-1]
#     obconv = openbabel.OBConversion()
#     obmol = openbabel.OBMol()
#     obconv.SetInFormat(ftype)
#     obconv.ReadFile(obmol,xyzf)
#     formula = obmol.GetFormula()
#     fout.write(formula + ',' + xyzf_)

# import sys,csv

# txt = sys.argv[1]
# f = open(txt,'r')
# csv

# ## ANN
# from molSimplify.Classes.mol3D import *
# from molSimplify.Classes.ligand import *
# from molSimplify.Scripts.tf_nn_prep import *
# import glob,csv

# xyz_path = '/Users/tzuhsiungyang/Dropbox (MIT)/Kulik Group Project Sean/all_oct_Fe_components_062618/'
# fpath = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/'
# fname = 'refcode_selected.txt'
# fout_name = 'eu_dist_selected.txt'
# f = open(fpath + fname,'r')
# for ref in f.readlines():
#     ref_ = ref.rstrip('\n')
#     xyzf = glob.glob(xyz_path + ref_ + '*xyz')[0]
#     xyz = mol3D()
#     xyz.readfromxyz(xyzf)
#     liglist, ligdents, ligcons = ligand_breakdown(xyz)
#     ax_ligand_list, eq_ligand_list, ax_natoms_list, eq_natoms_list, ax_con_int_list, eq_con_int_list, ax_con_list, eq_con_list, built_ligand_list = ligand_assign(xyz,liglist,ligdents,ligcons)
#     predictor = 'split'
#     # descriptor_names, descriptors = get_descriptor_vector(xyz)
#     metal = 'Fe'
#     oxidation_state = 2
#     result_dictionary = invoke_ANNs_from_mol3d(xyz,oxidation_state,alpha=0.2)
#     # result,latent_space_vector = ANN_supervisor(predictor,descriptors,descriptor_names)
#     # min_dist = find_true_min_eu_dist(predictor, descriptors, descriptor_names)
#     # print(min_dist)
#     dist = result_dictionary['distance']
#     dist.insert(0,ref_)
#     fout = open(fpath + fout_name,'a')
#     csv_out = csv.writer(fout)
#     csv_out.writerow(dist)
#     fout.close()

fin = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/Fe-N_BL.txt'



from scipy.optimize import curve_fit
from scipy.signal import argrelextrema

import matplotlib.pyplot as plt
import numpy as np

def hist_fitter(fin):

    # Make histograms

    f = open(fin,'r')
    a = []
    for line in f.readlines():
        line_ = line.rstrip('\n')
        a.append(float(line_))

    signal = False
    nbins = 50

    while not signal:
        a = np.array(a)
        plt.clf()
        n,bins,patches = plt.hist(a,nbins)
        x_range = [bins[0],bins[-1]+(bins[1]-bins[0])]
        bins = np.delete(bins,-1)
        y = n
        x = bins
        maxima = argrelextrema(y,np.greater)[0].tolist()
        if len(maxima) < 4:
            signal = True
        else:
            nbins -= 1
            print nbins

    # Gaussian fitting

    def func(x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            ctr = params[i]
            amp = params[i+1]
            wid = params[i+2]
            y = y + amp * np.exp(-((x - ctr)/wid)**2)
        return y

    guess = []
    for i in range(2):
        idx = maxima.pop()
        guess.append(x[idx])
        guess.append(y[idx])
        guess.append((x[idx+1]-x[idx-1])/1.5)

    lower_bound = [i*0.5 for i in guess]
    upper_bound = [i*4 for i in guess]

    popt, pcov = curve_fit(func, x, y, p0=guess,bounds=(lower_bound,upper_bound))
    print popt

    eqn = 'y = a*exp(-1/2*((x-b)/c)^2)'
    first_g = 'a = ' + str(round(popt[0],2)) + '\n'
    first_g += 'b = ' + str(round(popt[1],0)) + '\n'
    first_g += 'c = ' + str(round(popt[2],3))
    second_g = 'a = ' + str(round(popt[3],2)) + '\n'
    second_g += 'b = ' + str(round(popt[4],0)) + '\n'
    second_g += 'c = ' + str(round(popt[5],3))

    # # Mn-N text box
    # plt.text(1.9,800,eqn,bbox=dict(alpha=0.5))
    # plt.text(2.0,380,second_g)
    # plt.text(2.15,600,first_g)

    fontsize = 16

    # Fe-N textbox
    plt.text(1.65,1200,eqn,bbox=dict(alpha=0.5),fontsize=fontsize,fontweight='bold')
    plt.text(1.80,900,second_g,fontsize=fontsize,fontweight='bold')
    plt.text(2.12,600,first_g,fontsize=fontsize,fontweight='bold')

    # # Co-N textbox
    # plt.text(2.1,2800,eqn,bbox=dict(alpha=0.5))
    # plt.text(1.83,2200,second_g)
    # plt.text(2.12,1300,first_g)

    font = {'size':'22','weight':'bold','xtick':'22','ytick':'22'}
    xlabel_ = 'Fe-N distance / A'
    ylabel_ = 'Frequency'
    # plt.plot(x, y)
    x_smooth = np.linspace(x[0],x[-1],num=100)
    fit = func(x_smooth, *popt)
    # plt.rc('font',**font)
    plt.xlim(x_range)
    plt.plot(x_smooth, fit , 'r-')
    plt.xticks(fontsize=fontsize,fontweight='bold')
    plt.yticks(fontsize=fontsize,fontweight='bold')
    plt.xlabel(xlabel_,fontsize=fontsize,fontweight='bold')
    plt.ylabel(ylabel_,fontsize=fontsize,fontweight='bold')
    plt.show()
    # x_lower_upper = popt[0]-popt[2]*3
    # x_fill = np.linspace(x[0],x_lower_upper,num=100)
    # y_fill = func(x_fill, *popt)
    # plt.fill_between(x_fill,y_fill,0,facecolor='r')
    # x_fill_list.pop()
    # x_upper_lower = popt[3]+popt[5]*3
    # x_fill2 = np.linspace(x_upper_lower,x[-1],num=100)
    # y_fill2 = func(x_fill2, *popt)
    # plt.fill_between(x_fill2,y_fill2,0,facecolor='r')

hist_fitter(fin)

## Keras

# Generate dummy data
from keras.models import Sequential
from keras.layers import Dense, Activation
import numpy as np
fin = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/ref_for_sco_xyz_features.txt'

with open(fin,'r') as f1:
    d = []
    l = []
    i = 0
    for line in f1.readlines():
        if i > 0:
            line_ = line.rstrip('\n\r')
            d.append([float(j) for j in line_.split(',')[2:]])
            if str(line_.split(',')[1]) == 'hs':
                l.append(int(1))
            else:
                l.append(int(0))
        i += 1

data = np.asarray(d)
labels = np.asarray(l)

# For a single-input model with 2 classes (binary classification):
dim = data.shape[1]
model = Sequential()
model.add(Dense(32, activation='relu', input_dim=dim))
model.add(Dense(1, activation='sigmoid'))
# model.compile(optimizer='rmsprop',
#     loss='binary_crossentropy',
#     metrics=['accuracy'])
model.compile(optimizer='rmsprop',
    loss='binary_crossentropy',
    metrics=['accuracy'])

# Train the model, iterating on the data in batches of 32 samples
model.fit(data, labels, epochs=1000, batch_size=32)

## scikit learn SGDClassifier
import numpy as np
from sklearn import linear_model
from sklearn.preprocessing import normalize
# X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
# Y = np.array([1, 1, 2, 2])

feature_depth = 33
x = []
y = []
z = []
m = []
ftrain = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_features_with_label.txt'
f1 = open(ftrain,'r')
# for i,line in enumerate(f1.readlines()):
#     line_ = line.rstrip('\r\n')
#     line_ls = line_.split(',')
#     if i == 0:
#         header = line_ls[3:]
#     if i > 0:
#         atno = int(line_ls[2])
#         if atno == 26:
#             m.append(line[3:feature_depth])

for i,line in enumerate(f1.readlines()):
    line_ = line.rstrip('\r\n')
    line_ls = line_.split(',')
    if i == 0:
        header = line_ls[3:]
    if i > 0:
        atno = int(line_ls[2])
        if atno == 26:
            x.append([float(i) for i in line_ls[3:feature_depth]])
            if str(line_ls[1]) == 'hs':
                y.append(int(1))
            else:
                y.append(int(0))
            z.append(line_ls[0])

data = normalize(np.asarray(x).T).T
labels = np.asarray(y)

clf = linear_model.SGDClassifier(max_iter=1000)
clf.fit(data, labels)

x = []
y = []
z = []
ftest = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/ref_for_sco_xyz_features.txt'
f2 = open(ftest,'r')
for i,line in enumerate(f2.readlines()):
    line_ = line.rstrip('\r\n')
    line_ls = line_.split(',')
    if i == 0:
        header = line_ls[2:]
    if i > 0:
        atno = int(line_ls[2])
        if atno > 21:
            x.append([float(i) for i in line_ls[2:feature_depth]])
            if str(line_ls[1]) == 'hs':
                y.append(int(1))
            else:
                y.append(int(0))
            z.append(line_ls[0])

x_test = np.asarray(x)
y_test = np.asarray(y)
y_predict = clf.predict(x_test)
# fin2 = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_features_with_label.txt'
# f2 = open(fin2,'r')
# print(f2.readlines()[3044].rstrip('\r\n').split(',')[0])

## scikit-learn decision tree
import numpy as np
from sklearn import tree
from sklearn.preprocessing import normalize
import graphviz
# X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
# Y = np.array([1, 1, 2, 2])
# Data
feature_depth = 24
x = []
y = []
z = []
fin2 = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_features_with_label.txt'
f2 = open(fin2,'r')
for i,line in enumerate(f2.readlines()):
    line_ = line.rstrip('\r\n')
    line_ls = line_.split(',')
    if i == 0:
        header = line_ls[9:]
    if i > 0:
        atno = int(line_ls[2])
        if atno == 26:
            x.append([float(i) for i in line_ls[9:feature_depth]])
            if str(line_ls[1]) == 'hs':
                y.append(int(1))
            else:
                y.append(int(0))
            z.append(line_ls[0])

data = normalize(np.asarray(x).T).T
labels = np.asarray(y)
# Training
clf = tree.DecisionTreeClassifier(criterion='entropy',max_depth=5)
clf.fit(data, labels)
# Graph
dot_data = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dot_data)
graph.render('/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/iris')

## Tree-based feature selection
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import tree
from sklearn.feature_selection import SelectFromModel
import pandas as pd
import graphviz
# X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
# Y = np.array([1, 1, 2, 2])
# Data
feature_depth = 27
x = []
y = []
z = []
fin2 = '/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/all_oct_1TM_features_with_label.txt'
data_ = pd.read_csv(fin2)
data_ = data_.iloc[:,0:9]
data_ = pd.get_dummies(data_,columns=['mato','fatno_1','fatno_2','fatno_3','fatno_4','fatno_5','fatno_6'])
#,'fval_1','fval_2','fval_3','fval_4','fval_5','fval_6'])
headers = data_.columns[2:].values
labels = data_['ss'].values
data = data_.iloc[:,2:].values

# f2 = open(fin2,'r')
# for i,line in enumerate(f2.readlines()):
#     line_ = line.rstrip('\r\n')
#     line_ls = line_.split(',')
#     if i == 0:
#         header = line_ls[2:feature_depth]
#     if i > 0:
#         atno = int(line_ls[2])
#         if atno > 21:
#             x.append([float(i) for i in line_ls[2:feature_depth]])
#             if str(line_ls[1]) == 'hs':
#                 y.append(int(1))
#             else:
#                 y.append(int(0))
#             z.append(line_ls[0])

# data = normalize(np.asarray(x).T).T
# data = np.asarray(x)
# labels = np.asarray(y)
# data.shape

# clf = ExtraTreesClassifier(n_estimators=50)
# clf = clf.fit(data,labels)
# clf.feature_importances_
# model = SelectFromModel(clf,prefit=True)
# data_new = model.transform(data)
# data_new.shape

clf = tree.DecisionTreeClassifier(criterion='entropy',max_depth=5)
clf.fit(data, labels)
# Graph
dot_data = tree.export_graphviz(clf, out_file=None)
graph = graphviz.Source(dot_data)
graph.render('/Users/tzuhsiungyang/Dropbox (MIT)/Work at the Kulik group/bond length project/iris')

from astropy.table import Table,Column

t = Table( )