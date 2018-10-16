## extract hessian from .hess and then diagonalize

import glob, os, re
import numpy as np
# from molSimplify.Classes.mol3D import *

# partial hessian

partial_hess = [[0.041677998743, 0.041179646923, 0.011033922837, 0.0012780504069, 0.0023894512887, 0.0043170271304, -0.026985523574, -0.043959595491, -0.015853455845],
    [0.041179646923, 0.1116341609, 0.0078639300402, 0.0059849224541, -0.015229505848, 0.004878834599, -0.052726610788, -0.074365135965, -0.029028252688], 
    [0.011033922837, 0.0078639300402, 0.09599160515, 0.027114763354, 0.027599925377, 0.00099830785582, -0.041304880141, -0.054636509255, -0.017319539258], 
    [0.0012780504069, 0.0059849224541, 0.027114763354, -0.0069679064636, -0.024128860561, -0.026611331339, 0.0152634557, 0.033455851278, 0.01940617514], 
    [0.0023894512887, -0.015229505848, 0.027599925377, -0.024128860561, -0.010867956478, -0.033130689958, 0.045632795266, 0.051813268426, 0.032639257168], 
    [0.0043170271304, 0.004878834599, 0.00099830785582, -0.026611331339, -0.033130689958, 0.018384216132, 0.037284187187, 0.0503594692, 0.0010729231416], 
    [-0.026985523574, -0.052726610788, -0.041304880141, 0.0152634557, 0.045632795266, 0.037284187187, 0.43100707883, -0.19705412973, -0.072901745246], 
    [-0.043959595491, -0.074365135965, -0.054636509255, 0.033455851278, 0.051813268426, 0.0503594692, -0.19705412973, 0.27257520281, -0.095167291695], 
    [-0.015853455845, -0.029028252688, -0.017319539258, 0.01940617514, 0.032639257168, 0.0010729231416, -0.072901745246, -0.095167291695, 0.51996030059]]

# extract hessian from .hess
hess = glob.glob('*hess')[0]
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
    for j in hessian[1+i*(dof+1):(dof+1)+i*(dof+1)]:
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
idx_list = [27,28,40]
dof = [range(idx*3,idx*3+3) for idx in idx_list]
a = 0
for i in range(len(dof)):
    for ii in range(len(dof[i])):
        b = 0
        idx_i = dof[i][ii]
        for j in range(len(dof)):
            for jj in range(len(dof[j])):
                idx_j = dof[j][jj]
                npf[idx_i,idx_j] = partial_hess[a][b]
                b += 1
        a += 1

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