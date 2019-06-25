# -*- coding: utf-8 -*-
"""
Created on Fri May 25 16:11:03 2018

@author: kncv078
"""
# NOTES:
# fingerprint files are very large, i should use sparse matrices


# 1 - import libraries

import os
import time
import math
#import random
import pandas as pd
#import numpy as np
from scipy import sparse
import matplotlib as mpl
mpl.use('Agg') # this allows plt functions to work on scp even tough no display possible
from sklearn.ensemble import RandomForestClassifier
#from sklearn.metrics import roc_auc_score, matthews_corrcoef, make_scorer
from sklearn import metrics
import matplotlib.pyplot as plt
#from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import StratifiedKFold
#import my_functions as mf
#from sklearn.model_selection import train_test_split
#import seaborn
#from sklearn.svm import SVC

def get_metrics(labels_all, pred_all, Pred_prob_A_all,compounds):
    (tp,fn),(fp,tn) = metrics.confusion_matrix(labels_all, pred_all)
    A_count = tp+fn
    N_count = tn+fp
    mcc = ((tp*tn)-(fp*fn))/math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    F1 = (2*tp)/(2*tp+fp+fn)
    kappa = metrics.cohen_kappa_score(labels_all,pred_all)
    print('\n\n#_cmpds:',len(compounds),'\ttotal A:',A_count,'\ttotal N:',N_count)
    print('True Positive:\t', tp)
    print('False Negative:\t', fn)
    print('False Positive:\t', fp)
    print('True Negative:\t', tn)
    print('MCC:\t\t{}'.format(mcc))
    print('Kappa:\t\t{}'.format(kappa))
    if (tp or fn or fp) > 0:
        print('Precision:\t', precision)
        print('Recall:\t\t', recall)
        print('F1 score:\t', F1)
    fpr, tpr, thresholds = metrics.roc_curve(labels_all, Pred_prob_A_all, pos_label='A')
    roc_auc = metrics.auc(fpr, tpr)
    print('ROC AUC:\t', roc_auc)
    return (tp, fn, tn, fp, mcc, precision, recall, F1, kappa, fpr, tpr, roc_auc)

def plot_roc(fpr, tpr, roc_auc, outpath, assay, name, cvf):
    plt.figure()
    plt.plot(fpr, tpr, lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1],[0, 1], lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC - {}-{}'.format(assay,name))
    plt.legend(loc="lower right")
    plt.savefig('{}{}_{}_{}foldCV_ROC_curve.png'.format(outpath,assay,name,cvf), bbox_inches='tight')
    #plt.show()
    plt.close()

def get_assay_list(assay_list_file,ntests):
    print('-loading assay list')
    assay_list =[]
    with open(assay_list_file, 'r') as inf:
        for line in inf:
            assay_list.append(str(line.strip()))
    #test_assays=random.sample(assay_list,ntests)
    test_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553',
                   '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
#    test_assays = ['798','1515','2553','463104','504454','588497','624414','686964'] # original 8 validation assays
#    test_assays = ['2732','1236','720543','588621','463165','602229','1899','834','1510','560'] #for HP search
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
    print('\nloading FP df...')
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

def main1_():
    start = time.time()
    outdir = "CrossValidation2/Assay_{}/"
    ntests = 20 #number of test assays to investigate - randomly chosen
    CVfolds = 6 #number of CV folds to test
    n_threads = 6
    nTrees = 150
    Ranstate = 56


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
    else: print('WARNING: compound lists do not match!')

    # In[] starting main part

    print('df size:', len(HTSFP))
    FP_dict = {'htsfp':binhtsfp_file,'ecfp':ecfp4_file,'cesfp':(binhtsfp_file,ecfp4_file)}
    for FPname, FPpath in FP_dict.items():
#        if FPname == 'ecfp': continue
        FP = get_sparse_fp(FPpath,FPname)
        print('###'*12)
        print('### Beginning analysis of: {} ###'.format(FPname))
        print('###'*12)
        for assay in test_assaysD.keys():
            print('---'*12)
            print('Assay:',assay)
            outpath = outdir.format(assay)

            assay_idx = assayls.index(assay)
            label = HTSFP[assay_idx][HTSFP[assay_idx] !='x']

            # 14 - check to see if selected new assay has all flags
            skip, As = get_xNA_counts(label,CVfolds)
            if skip: continue

            req_idxs = set(str(c) for c in label.index)
            cmpd_idx_list = [i for i,j in enumerate(cmpdls) if j in req_idxs]

            #     cmpd_idx_list.append(cmpdls.index(str(cmpd)))
            # 10 - determine intersection of compound lists newAssay_cmpds and new11Assays
            print('number of compounds in dataset:\t\t{}'.format(FP.shape[0]))
            print('number of compounds in assay {}:\t{}'.format(assay, len(label)))

            ac_fp = FP[cmpd_idx_list]
            if FPname == 'htsfp':
                if int(ac_fp[:,assay_idx].count_nonzero()) != As: 
                    print('labels dont match extracted sparse column',int(ac_fp[:,assay_idx].count_nonzero()),As)
                keep_cols = list(range(561))
                del keep_cols[assay_idx]
                ac_fp = ac_fp[:,keep_cols]
                rf = RandomForestClassifier(n_jobs=n_threads, n_estimators=nTrees, class_weight='balanced',
                                            max_features='sqrt', criterion='entropy', max_depth=40, min_samples_split=2,
                                            min_samples_leaf=5, random_state=Ranstate)
            elif FPname == 'cesfp':
                if int(ac_fp[:,assay_idx].count_nonzero()) != As: 
                    print('labels dont match extracted sparse column',int(ac_fp[:,assay_idx].count_nonzero()),As)
                keep_cols = list(range(561+1024))
                del keep_cols[assay_idx]
                ac_fp = ac_fp[:,keep_cols]
                rf = RandomForestClassifier(n_jobs=n_threads, n_estimators=nTrees, class_weight='balanced',
                                            max_features='sqrt', criterion='gini', max_depth=None, min_samples_split=2,
                                            min_samples_leaf=8, random_state=Ranstate)
            elif FPname == 'ecfp':
                rf = RandomForestClassifier(n_jobs=n_threads, n_estimators=nTrees+50, class_weight='balanced',
                                            max_features='sqrt', criterion='gini', max_depth=30, min_samples_split=2,
                                            min_samples_leaf=8, random_state=Ranstate)

            ## RANDOM FOREST ANALYSIS
            ## results from grid search:
#            ecfp HP settings
#            'n_estimators': 150, 'class_weight': 'balanced', 'oob_score': True, 'max_features':'sqrt', 'criterion': 'gini', 'max_depth': 30, 'min_samples_split': 8, 'min_samples_leaf': 3
#            htsfp HP settings
#            'n_estimators': 150, 'class_weight': 'balanced', 'oob_score': True, 'max_features':'sqrt', 'criterion': 'entropy', 'max_depth': 20, 'min_samples_split': 6, 'min_samples_leaf': 3
#            BaSH HP settings
#            'n_estimators': 150, 'class_weight': 'balanced', 'oob_score': True, 'max_features':'sqrt', 'criterion': 'entropy', 'max_depth': 20, 'min_samples_split': 8, 'min_samples_leaf': 3
            compounds = []
            labels_all = []
            Pred_prob_A_all = []
            pred_all = []
            
            if os.path.exists('{}{}_{}_cvfold{}_predictions.txt'.format(outpath,assay,FPname,CVfolds)) == True:
                print('file with name: {} already exists under this outpath, moving on to next assay'.format(assay))
                continue

            print('Starting random forest')
            skf = StratifiedKFold(n_splits=CVfolds, random_state=Ranstate, shuffle=False)
            cvfld = 0
            for train_index, test_index in skf.split(ac_fp, label):
                fp_train, fp_test = ac_fp[train_index], ac_fp[test_index]
                L_train, L_test = label.iloc[train_index], label.iloc[test_index]
                cvfld+=1
                print('Fold {}...  '.format(cvfld),end='\r')

#                print('fitting...', end='\t')
                rf.fit(fp_train,L_train)
#                print('predicting...', end='\t')
                pred = rf.predict(fp_test)
#                print('predicting with probability scores...')
                pred_prob = rf.predict_proba(fp_test)
                Feat_imp = rf.feature_importances_
            
                if FPname == 'cesfp':
                    with open('FI/Feature_importance_{}_{}_{}.out'.format(FPname,assay,cvfld), 'w') as f:
                        for i in Feat_imp:
                            f.write('{}\n'.format(i))

                compounds += list(L_test.index)
                labels_all += list(L_test)
                Pred_prob_A_all += list(pred_prob[:,0])
                pred_all += list(pred)
                
                pp = pd.DataFrame({'cmpd': list(L_test.index),'label':list(L_test),'pred':list(pred), 'PredProb(A)': list(pred_prob[:,0])})
                pp = pp.set_index('cmpd')
                pp = pp.sort_values('PredProb(A)', ascending=False)
                pp.to_csv('{}{}_{}_cvfold{}_predictions.txt'.format(outpath,assay,FPname,cvfld), sep='\t')
            #save files to outpath directory
            #create new dir if not already existing
            if not os.path.exists(outpath):
                os.makedirs(outpath)

            ppf = pd.DataFrame({'cmpd': compounds,'label':labels_all,'pred':pred_all, 'PredProb(A)': Pred_prob_A_all})
            ppf = ppf.set_index('cmpd')
            ppf = ppf.sort_values('PredProb(A)', ascending=False)
            ppf.to_csv('{}{}_{}_full_predictions.txt'.format(outpath,assay,FPname), sep='\t')
            
            tp, fn, tn, fp, mcc, precision, recall, F1, kappa, fpr, tpr, roc_auc = get_metrics(labels_all, pred_all, Pred_prob_A_all, compounds)
            plot_roc(fpr, tpr, roc_auc, outpath, assay, FPname, cvfld)

    print('timetaken: {}'.format(time.time()-start))
main1_()


