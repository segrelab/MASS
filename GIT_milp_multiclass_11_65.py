# -*- coding: utf-8 -*-
"""Prediction of representative phenotypes using multi-output subset selection
"""

from gurobipy import *
import time
import numpy as np
import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from numpy import linalg as LA

from sklearn.metrics import confusion_matrix,classification_report,f1_score,accuracy_score,r2_score
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.model_selection import cross_val_score
from numpy.linalg import inv
from iterstrat.ml_stratifiers import MultilabelStratifiedKFold
import warnings
warnings.filterwarnings("ignore")


from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import LassoCV
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
# scaler = StandardScaler()
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier



def MMRP_classify2021(X, r, z_start=None, z_greedy=None,n_greedy=None, OutputFlag=0,Ncpu=30,LP=False,lamda=0,include=None,exclude=None,MIPFocus=0):#,M=100
    n, d = X.shape   
    # assert n == len(y)
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
        
    m.addConstr(quicksum(z[j] for j in range(d))<=r )
    for i in range(n):     
        for j in range(d): 
              # overloaded form
            m.addConstr((z[j] == 0) >> (X[i,j]*(b0[j]+quicksum([B[k,j]*X[i,k] for k in range(d)]))>=1-KSI[i,j]))


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
    # m.Params.method = 2
    #     m.params.DisplayInterval = 300
    m.Params.Threads = Ncpu
    m.Params.OutputFlag = OutputFlag#1#0#
    if MIPFocus!=0:
        m.Params.MIPFocus=MIPFocus
    #     m.setParam("TimeLimit", TimeLimit)
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
    
    return bestB, bestb0, bestz, bestKSI, m.getObjective().getValue()

def growth_categories(df):
  # cast symbols into numbers
  # - -> 0
    df.replace(to_replace='-',
             value=0.0,
             inplace=True)
    # V -> 1
    df.replace(to_replace='V',
             value=1.0,
             inplace=True)
    # W -> 1
    df.replace(to_replace='W',
             value=1.0,
             inplace=True)
    # D -> 1
    df.replace(to_replace='D',
             value=1.0,
             inplace=True)
    # + -> 2
    df.replace(to_replace='+',
             value=2.0,
             inplace=True)
    # ? -> NaN
    df.replace(to_replace='?',
             value=np.NaN,
             inplace=True)
    return df

def get_score(df_z_all, f1_micro_all, names_y_all, col_score='f1_score',clf='classify_RF'):
    output=f'../result/MIP_{clf}_{col_score}.csv'
    f1_scores=[]
    for i in range(len(f1_micro_all)):
        for nam, j in zip(names_y_all[i],f1_micro_all[i]):
            f1_scores.append((df_z_all.shape[0]-i,nam, j))
    f1_scores=pd.DataFrame.from_records(f1_scores, columns=['# of predictors','response',col_score])
    f1_scores.to_csv(output,index=False)
    return f1_scores

def plot_score(f1_scores, col_score='f1_score'):
    a4_dims = (14.7, 8.27)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.scatterplot(ax=ax, x="# of predictors", y=col_score, data=f1_scores)

    fig, ax = plt.subplots(figsize=a4_dims)
    sns.boxplot(ax=ax, x="# of predictors", y=col_score, data=f1_scores)

    fig, ax = plt.subplots(figsize=a4_dims)
    sns.pointplot(ax=ax, x="# of predictors", y=col_score, data=f1_scores,ci="sd")

if __name__ == '__main__':
    # M = 1000
    # Ncpu = 16
    # TimeLimit = 1500
    # Snitkin_RefinedExperimentalData.csv
    df=pd.read_csv('../data/growth_profiles_binary.csv',index_col='strain')  #,index_col='Mutant'

    df_new=df.astype(float)
    X=df_new.values
    n, d = X.shape
    """## only for classify, 0 -->-1 """

    X[X==0]=-1

    """# MMRP_classify

    ['D-Glucose_1.0', 'D-Glucose_2.0', 'Inulin_1.0', 'Inulin_2.0',
            'myo-Inositol_1.0', 'myo-Inositol_2.0', 'Methanol_1.0',
            'Methanol_2.0']
    """

    my_range=range(d-1,0,-1)#range(d-2,0,-2)#
    my_l=len(my_range)#d-1#37
    myseed = 2020
    rand = np.random.RandomState([myseed])

    z_all=[]
    B_all=[]
    obj_all=[]
    for r in my_range: 
        bestB, bestb0, bestz, bestKSI, obj = MMRP_classify2021(X, r)
        z_all.append(bestz)
        B_all.append(bestB)
        obj_all.append(obj)

    
    plt.figure(figsize=(20,8))
    sns.heatmap(np.array(z_all).astype(int), 
                xticklabels=df_new.columns.values,
                yticklabels=my_range, vmin=0, vmax=1)

    df_z_all=pd.DataFrame.from_records(z_all, columns=df_new.columns)
    df_z_all.to_csv('../result/MIP_classify_z_all_11.csv',index=False)
    """## obj_all increase when we have fewer features and more responses"""

    pd.Series(obj_all).plot()

    plt.figure(figsize=(20,18))
    sns.heatmap(np.array(z_all).T.astype(int), 
                xticklabels=my_range,#str_times,
                yticklabels=df_new.columns.values, vmin=0, vmax=1)

    df3sum=df_z_all.sum()
    df3sum.plot()

    """# times of selected"""

    df_z_all=pd.read_csv('../result/MIP_classify_z_all_11.csv')

    dfsum=df.sum().to_frame('sum of variables')
    dfsum['sum of selected times']=df_z_all.sum().values
    dfsum.sort_values('sum of variables').plot(x='sum of variables',y='sum of selected times')

    a4_dims = (6.7, 4.27)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.scatterplot(ax=ax, x="sum of variables", y="sum of selected times", data=dfsum.sort_values('sum of variables'))

    dfstd=df.std().to_frame('std of variables')
    dfstd['sum of selected times']=df_z_all.sum().values
    a4_dims = (6.7, 4.27)
    fig, ax = plt.subplots(figsize=a4_dims)
    sns.scatterplot(ax=ax, x="std of variables", y="sum of selected times", data=dfstd.sort_values('std of variables'))

    """# predict setting--MultilabelStratifiedKFold

    iterative-stratification is a project that provides scikit-learn compatible cross validators with stratification for multilabel data.

    Presently scikit-learn provides several cross validators with stratification. However, these cross validators do not offer the ability to stratify multilabel data. This iterative-stratification project offers implementations of MultilabelStratifiedKFold, MultilabelRepeatedStratifiedKFold, and MultilabelStratifiedShuffleSplit with a base algorithm for stratifying multilabel data described in the following paper:

    Sechidis K., Tsoumakas G., Vlahavas I. (2011) On the Stratification of Multi-Label Data. In: Gunopulos D., Hofmann T., Malerba D., Vazirgiannis M. (eds) Machine Learning and Knowledge Discovery in Databases. ECML PKDD 2011. Lecture Notes in Computer Science, vol 6913. Springer, Berlin, Heidelberg.
    """

    # dfc=pd.read_csv('../data/growth_profiles_continuous.csv',index_col='strain')  #,index_col='Mutant'

    # from sklearn.preprocessing import KBinsDiscretizer
    # est = KBinsDiscretizer(n_bins=10, encode='ordinal', strategy='quantile')
    # y_=est.fit_transform(dfc['HMBoligo'].values.reshape(-1, 1))
    # pd.Series(np.array(y_).ravel()).value_counts()

    X=df.values
    names_raw=df.columns.values
    n, d = X.shape  

    z_last = df_z_all.iloc[-1,:].values
    y_last = X[:,z_last==0] 


    myseed = 2020
    cvFolds= 5
    r_ = np.random.RandomState([myseed])
    """# LogisticRegressionCV"""

    acc_all=[]
    f1_micro_all=[]
    for i in range(0,df_z_all.shape[0]):
        z=df_z_all.iloc[i,:].values
        y_split = X[:,z==0]
        if i==0:
            folds = StratifiedKFold(n_splits=cvFolds, shuffle=True, random_state=myseed)
        else:
            folds = MultilabelStratifiedKFold(n_splits=cvFolds, shuffle=True, random_state=myseed)
        name_i = names_raw[z==0]
        for fold_, (trn_idx, test_idx) in enumerate(folds.split(X,y_split)):        
            print('r, and fold_==', d-i-1,fold_, end=';  ')  
            train_X, test_X = X[trn_idx][:,z==1], X[test_idx][:,z==1]
            train_y, test_y = X[trn_idx][:,z==0], X[test_idx][:,z==0]
            
            if fold_==0:
                print('name_i,train_y.shape, train_y.mean, test_y.mean:',name_i,train_y.shape, train_y.mean(), test_y.mean())
    #             train_X=scaler.fit_transform(train_X)
    #             test_X=scaler.transform(test_X)
            for j in range(train_y.shape[1]):     #penalty='l1', multi_class='auto',
                if (train_y[:,j]==0).sum()>cvFolds:
                    clf = LogisticRegressionCV(Cs=[10**-2,10**-1,10**0],class_weight='balanced', tol=0.001, cv=cvFolds, random_state=myseed).fit(train_X, train_y[:,j])
                else: # not cv due to ValueError: Class label 0 not present.
                    clf = LogisticRegression(C=10**-1,class_weight='balanced', tol=0.001,random_state=myseed).fit(train_X, train_y[:,j])

                y_te_pre=clf.predict(test_X)
                acc_all.append(accuracy_score(test_y[:,j],y_te_pre))
                f1_micro_all.append([df_z_all.shape[0]-i,fold_, name_i[j], f1_score(test_y[:,j],y_te_pre, average='micro')])

    f1_scores = get_score(f1_micro_all, col_score='f1_score',clf='classify_LR')
    f1_scores.tail()

    """# RandomForestClassifier"""

    acc_all=[]
    f1_micro_all=[]

    for i in range(0,df_z_all.shape[0]):
        z=df_z_all.iloc[i,:].values
        y_split = X[:,z==0]
        if i==0:
            folds = StratifiedKFold(n_splits=cvFolds, shuffle=True, random_state=myseed)
        else:
            folds = MultilabelStratifiedKFold(n_splits=cvFolds, shuffle=True, random_state=myseed)
        name_i = names_raw[z==0]

        for fold_, (trn_idx, test_idx) in enumerate(folds.split(X,y_split)):        
            print('r, and fold_==', d-i-1,fold_, end=';  ')  
            train_X, test_X = X[trn_idx][:,z==1], X[test_idx][:,z==1]
            train_y, test_y = X[trn_idx][:,z==0], X[test_idx][:,z==0]
            
            if fold_==0:
                print('name_i,train_y.shape, train_y.mean, test_y.mean:',name_i,train_y.shape, train_y.mean(), test_y.mean())
    #             train_X=scaler.fit_transform(train_X)
    #             test_X=scaler.transform(test_X)
            rfc=RandomForestClassifier(n_estimators=10, random_state=myseed,n_jobs=5, min_samples_split = 4, class_weight="balanced")

            for j in range(train_y.shape[1]):     #penalty='l1', multi_class='auto',
                if (train_y[:,j]==0).sum()>cvFolds:
                    rfc.fit(train_X, train_y[:,j])
                    y_te_pre=rfc.predict(test_X)
                else: # not cv due to ValueError: Class label 0 not present.
                    CV_rfc = GridSearchCV(estimator=rfc, param_grid=param_grid, cv= cvFolds,n_jobs=9)#, return_train_score=True iid=False,
                    CV_rfc.fit(train_X, train_y[:,j])
                    y_te_pre=CV_rfc.predict(test_X)
                acc_all.append(accuracy_score(test_y[:,j],y_te_pre))
                f1_micro_all.append([df_z_all.shape[0]-i,fold_, name_i[j], f1_score(test_y[:,j],y_te_pre, average='micro')])

    f1_scores = get_score(f1_micro_all, col_score='f1_score',clf='classify_RF')
    f1_scores.tail()

    """# compare with enumerate when 1 response"""
    acc_all=[]
    f1_micro_all=[]

    for i in range(0, X.shape[1]):
        X_split = X[:,[i_ for i_ in range(X.shape[1]) if i_ !=i]]
        y_split = X[:,i].reshape(-1,1)
        folds = StratifiedKFold(n_splits=cvFolds, shuffle=True, random_state=myseed)
        name_i = names_raw[i]

        for fold_, (trn_idx, test_idx) in enumerate(folds.split(X,y_split)):        
            print('r, and fold_==', d-i-1,fold_, end=';  ')  
            train_X, test_X = X_split[trn_idx], X_split[test_idx]
            train_y, test_y = y_split[trn_idx], y_split[test_idx]
            
            if fold_==0:
                print('name_i,train_y.shape, train_y.mean, test_y.mean:',name_i,train_y.shape, train_y.mean(), test_y.mean())

            rfc=RandomForestClassifier(n_estimators=10, random_state=myseed,n_jobs=5, min_samples_split = 4, class_weight="balanced")

            for j in range(train_y.shape[1]):     #penalty='l1', multi_class='auto',
                if (train_y[:,j]==0).sum()>cvFolds:
                    rfc.fit(train_X, train_y[:,j])
                    y_te_pre=rfc.predict(test_X)
                else: # not cv due to ValueError: Class label 0 not present.
                    CV_rfc = GridSearchCV(estimator=rfc, param_grid=param_grid, cv= cvFolds,n_jobs=9)#, return_train_score=True iid=False,
                    CV_rfc.fit(train_X, train_y[:,j])
                    y_te_pre=CV_rfc.predict(test_X)
                acc_all.append(accuracy_score(test_y[:,j],y_te_pre))
                f1_micro_all.append([i,fold_, name_i, f1_score(test_y[:,j],y_te_pre, average='micro')])

    f1_scores11 = get_score11(f1_micro_all, col_score='f1_score',clf='classify_RF')
    f1_scores11.plot(kind='bar',x='response', y='mean')

