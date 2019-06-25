# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:11:03 2018

@author: kncv078
"""
# NOTES:
# fingerprint files are very large, i should use sparse matrices


# 1 - import libraries

#import os
import time
import math
import random
import pandas as pd
#import numpy as np
from scipy import sparse
#import matplotlib as mpl
#mpl.use('Agg') # this allows plt functions to work on scp even tough no display possible
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, matthews_corrcoef, make_scorer
from sklearn import metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import GridSearchCV
#from sklearn.model_selection import StratifiedKFold
#import my_functions as mf
#from sklearn.model_selection import train_test_split
#import seaborn
#from sklearn.svm import SVC
    
def get_assay_list(assay_list_file,ntests):
    print('-loading assay list')
    assay_list =[]
    with open(assay_list_file, 'r') as inf:
        for line in inf:
            assay_list.append(str(line.strip()))
#    test_assays=random.sample(assay_list,ntests)
    test_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553',
                   '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700'] #validation
#    test_assays = ['2732','1236','720543','588621','463165','602229','1899','834','1510','560'] #for HP search
#    test_assays = ['504812','504810']
    print(test_assays)
    test_assay_dict = {}
    for a in test_assays:
        test_assay_dict[a] = assay_list.index(a)
    return test_assay_dict,assay_list

def get_labels(htsFP_file,testassays):
    print('-loading htsfp labels')
    st = time.time() 
    cols = [i+1 for i in testassays]
    print(cols)
    df = pd.read_csv(htsFP_file, sep=' ', header=None, index_col=0, usecols=[0]+cols, low_memory=True)
    df.columns = testassays
    print('time taken:',time.time()-st)
    return df

def get_xNA_counts(labs,folds):
    skip = False
    x_count = labs[labs=='x'].count()
    A_count = labs[labs=='A'].count()
    length = labs[labs!='x'].count()
    print('x_count:\t{}\nA_count:\t{}\nTotal(A&N):\t{}'.format(x_count, A_count, length))
    if A_count < folds: print('Warning: not enough actives in this assay'); skip=True
    return skip, A_count

def get_sparse_fp(path,name):
    print('loading FP df...')
    if name == 'cesfp':
        cooMecfp = sparse.load_npz(path[1])
        cooMhts = sparse.load_npz(path[0])
        cooM = sparse.hstack([cooMhts,cooMecfp])
    else:
        cooM = sparse.load_npz(path)
    return cooM.tocsr()

def load_cmpd_list(path):
    cmpdls = []
    with open(path,'r') as inf:
        for line in inf:
            cmpdls.append(line.strip())
    return cmpdls

def gridsearch(features,labels,folds,FPname,assay):    
    print('grid search...')
    parameters = {'n_estimators':[40],'criterion':['gini'],'class_weight':['balanced'],
                  'max_features':['sqrt'],'min_samples_leaf':[8], 'max_depth':[None], 
                  'min_samples_split':[2],'bootstrap':[True]}
    rf = RandomForestClassifier(n_jobs=6, random_state=56)
    clf = GridSearchCV(estimator=rf, param_grid=parameters, cv=folds, scoring=make_scorer(roc_auc_score, needs_proba=True))
    clf.fit(features, labels)
    print('{}:\t{}\nbestscore{}'.format(FPname,clf.best_params_,clf.best_score_))
    with open('C:/CESFP_project/HP_search.txt','a') as f:
        f.write('{}\t{}:\t{}\tbest score: {}\n'.format(assay,FPname,clf.best_params_,clf.best_score_))

def main1_():
    start = time.time()
    outdir = "CrossValidation/Assay_{}/"
    ntests = 561 #number of test assays to investigate
    CVfolds = 3 #number of CV folds to test
    
    
    # 2 - define input file paths
    assay_list_file= 'C:/CESFP_project/htsfp_sparse_assaylist.txt'
    htscmpdls_path = 'C:/CESFP_project/htsfp_sparse_cmpdlist.txt'
    ecfpcmpdls_path= 'C:/CESFP_project/ecfp_sparse_cmpdlist.txt'
    binhtsfp_file  = 'C:/CESFP_project/htsfp_sparse.npz'
    ecfp4_file     = 'C:/CESFP_project/ecfp_sparse.npz'
    htsFP_file     = 'C:/CESFP_project/step2_out_htsfp_t20000.txt'
    
    # In[1] - Load all required files.
    ecfpcmpdls = load_cmpd_list(ecfpcmpdls_path)
    htscmpdls = load_cmpd_list(htscmpdls_path)
    
    test_assaysD, assayls = get_assay_list(assay_list_file,ntests)
    HTSFP = get_labels(htsFP_file,sorted(test_assaysD.values()))
    
    match = True
    for i,j,k in zip(ecfpcmpdls,htscmpdls,HTSFP.index):
        if i != j or i!=str(k): match=False; break
    if match: 
        cmpdls=htscmpdls.copy()
        del ecfpcmpdls
        del htscmpdls
    else: print('WARNING: compund lists do not match!')
    
    # In[] starting main part
    
    print('df size:', len(HTSFP))
    FP_dict = {'htsfp':binhtsfp_file,'ecfp':ecfp4_file,'cesfp':(binhtsfp_file,ecfp4_file)}
    cz = []
    for FPname, FPpath in FP_dict.items():
        if FPname != 'cesfp': continue
        FP = get_sparse_fp(FPpath,FPname)
        print('###'*12)
        print('Beginning analysis of: {}'.format(FPname))
        for assay in test_assaysD.keys():
            print('---'*12)
            print('Assay:',assay)
            #outpath = outdir.format(assay)
            
            assay_idx = assayls.index(assay)
            label = HTSFP[assay_idx][HTSFP[assay_idx] !='x']
            
            # 14 - check to see if selected new assay has all flags
            skip, As = get_xNA_counts(label,CVfolds)
            if skip: continue
    
            req_idxs = set(str(c) for c in label.index)
            cmpd_idx_list = [i for i,j in enumerate(cmpdls) if j in req_idxs]
    
    #            cmpd_idx_list.append(cmpdls.index(str(cmpd)))
            # 10 - determine intersection of compound lists newAssay_cmpds and new11Assays
            print('number of compounds in dataset:\t\t{}'.format(FP.shape[0]))
            print('number of compounds in assay {}:\t{}'.format(assay, len(label)))
            
            
            ac_fp = FP[cmpd_idx_list]
            cz.append(int(ac_fp[:,assay_idx].count_nonzero())/len(label))
            if int(ac_fp[:,assay_idx].count_nonzero()) != As: print('labels dont match extracted sparse column',int(ac_fp[:,assay_idx].count_nonzero()),As)
            if FPname == 'htsfp':
                keep_cols = list(range(561))
                del keep_cols[assay_idx]
                ac_fp = ac_fp[:,keep_cols]
            elif FPname == 'cesfp':
                keep_cols = list(range(561+1024))
                del keep_cols[assay_idx]
                ac_fp = ac_fp[:,keep_cols]
            
            #ac_fp = ac_fp.drop([assay_idx],axis=1)
    
            ## RANDOM FOREST ANALYSIS
            with open('stuffiwant.csv','a') as outf:
                outf.write('{},{},{}\n'.format(assay,len(label),As))
            continue
            gridsearch(ac_fp, label, CVfolds,FPname,assay)
               
#    print('timetaken: {}'.format(time.time()-start))
main1_()


