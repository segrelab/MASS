# -*- coding: utf-8 -*-
"""
Prediction of representative phenotypes using multi-output subset selection
Author: wty@bu.edu; with minor clean-ups by Konrad Herbst, herbstk@bu.edu
If you find it useful, please cite our paper as follows:
@article{forchielli2022prediction,
  title={Prediction of representative phenotypes using multi-output subset selection},
  author={Forchielli, Elena and Wang, Taiyao and Thommes, Meghan and Paschalidis, Ioannis Ch and Segre, Daniel},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
@phdthesis{wang2020data,
  title={Data analytics and optimization methods in biomedical systems: from microbes to humans},
  author={Wang, Taiyao},
  year={2020}
}
"""

import numpy as np
import pandas as pd
from gurobipy import *

def MASS(X, r, z_start=None, z_greedy=None, n_greedy=None,
         OutputFlag=0, Ncpu=70, LP=False, lamda=0,
         include=None, exclude=None,
         MIPFocus=0, TimeLimit = None, pairs_ind=None):#,M=100
    n, d = X.shape
    m = Model()

    ## addVars
    B = m.addVars(d, d, lb = -GRB.INFINITY) 
    Babs = m.addVars(d, d, lb = 0) 
    b0 = m.addVars(d, lb = -GRB.INFINITY) #intercept
    z = m.addVars(d, vtype = GRB.BINARY)
    zz = m.addVars(d, vtype = GRB.BINARY)
#     Y = m.addVars(n, d, lb = -GRB.INFINITY) #t_ij
    KSI = m.addVars(n, d, lb = 0)
#     Y_X = m.addVars(n, d, lb = -GRB.INFINITY) 

    ## addConstr  
    # In an SOS constraint of type 1 (an SOS1 constraint), at most one variable in the specified list is allowed to take a non-zero value.
    for j in range(d):
        m.addConstr(zz[j]==1-z[j])
#     for j in range(0,d,2):# select group or not
#         m.addConstr(z[j]==z[j+1])
    if pairs_ind is not None:
        for pair in pairs_ind:
            p1,p2=pair
            m.addConstr(z[p1]==z[p2])
    m.addConstr(quicksum(z[j] for j in range(d))<=r )
    for i in range(n):     
        for j in range(d): 
              # overloaded form
            m.addConstr((z[j] == 0) >> (X[i,j]*(b0[j]+quicksum([B[k,j]*X[i,k] for k in range(d)]))>=1-KSI[i,j]))
#             m.addConstr(X[i,j]*(b0[j]+quicksum([B[k,j]*X[i,k] for k in range(d)]))>=1-KSI[i,j]-M*z[j])
#             m.addConstr(Y_X[i,j] == Y[i,j]-X[i,j])
#             m.addConstr( W[i,j]- (Y[i,j]-b0[j]-quicksum([B[k,j]*X[i,k] for k in range(d)])) >= 0)
#             m.addConstr( W[i,j]+ (Y[i,j]-b0[j]-quicksum([B[k,j]*X[i,k] for k in range(d)])) >= 0)
#             m.addSOS(GRB.SOS_TYPE1, [Y_X[i,j], zz[j]], [1, 2])

    for k in range(d): 
        for j in range(d): 
            m.addSOS(GRB.SOS_TYPE1, [B[j,k], zz[j]], [1, 2])
            m.addConstr(Babs[j,k] == abs_(B[j,k]), name="absconstr")
    if z_start is not None:
        for j in range(d):
            z[j].start = z_start[j]
            zz[j].start = 1-z_start[j]
    if z_greedy is not None:
        m.addConstr(quicksum(z_greedy[j]*z[j] for j in range(d)) >= r-n_greedy)
    if LP:
        for j in range(d):
            m.addConstr(z[j] == z_start[j])
            m.addConstr(zz[j]== 1-z_start[j])
    if include is not None:
        for j in include:
            m.addConstr(z[j] == 1)
    if exclude is not None:
        for j in exclude:
            m.addConstr(z[j] == 0)
    # # B[j,j]=1, if z_j=0;
    # for k in range(d):
    #     for j in range(k+1,d):

    # m.Params.method = 2
    #     m.params.DisplayInterval = 300
    m.Params.Threads = Ncpu
    m.Params.OutputFlag = OutputFlag#1#0#
    if MIPFocus!=0:
        m.Params.MIPFocus=MIPFocus
    if TimeLimit is not None:
        m.setParam("TimeLimit", TimeLimit)
    obj2 = quicksum(KSI[i,j] for i in range(n) for j in range(d))/n + lamda*quicksum(Babs[j,k] for k in range(d) for j in range(d))
    m.setObjective(obj2, GRB.MINIMIZE)
    ## solve model
    m.optimize() 

    bestB = np.zeros((d,d))
    for i in range(d):
        for j in range(d):
            bestB[i, j] = B[i,j].X
    bestKSI = np.zeros((n,d))
    for i in range(n):
        for j in range(d):
            bestKSI[i, j] = KSI[i,j].X
    bestb0 = np.zeros((d,))
    bestz = np.zeros((d,))
    for j in range(d):
        bestb0[j] = b0[j].X
        bestz[j] = z[j].X

    return bestB, bestb0, bestz, bestKSI, m.getObjective().getValue(), m

def get_score(df_z_all, f1_micro_all, names_y_all, col_score='f1_score', clf='classify_RF'):
    output = f'../Results/MIP_{clf}_{col_score}.csv'
    f1_scores = []
    for i in range(len(f1_micro_all)):
        for nam, j in zip(names_y_all[i], f1_micro_all[i]):
            f1_scores.append((df_z_all.shape[0]-i, nam, j))
    f1_scores = pd.DataFrame.from_records(
        f1_scores, columns=['# of predictors', 'response', col_score])
    f1_scores.to_csv(output, index=False)
    return f1_scores

