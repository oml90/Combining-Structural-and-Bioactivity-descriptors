# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:07:08 2019

@author: laufkoetter
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 15:34:50 2018

@author: kncv078
"""

def do_maths(assay, fp, m):
    df = pd.read_csv(file_path.format(assay,assay,fp), sep='\t').sort_values('PredProb(A)', ascending=False).reset_index(drop=True)
    
    fpr1, tpr1, thresholds1 = metrics.roc_curve(df['label'], df['PredProb(A)'], pos_label='A')
    roc_auc1 = metrics.auc(fpr1, tpr1)
    roc_auc1 = round(roc_auc1,3)
    prf1 = metrics.precision_recall_fscore_support(df['label'], df['pred'])
    CM1 = metrics.confusion_matrix(df['label'], df['pred'])
    kappa1 = metrics.cohen_kappa_score(df['label'], df['pred'])
    MCC = metrics.matthews_corrcoef(df['label'], df['pred'])
    A_count1 = prf1[3][0]
    N_count1 = prf1[3][1]
    P1 = prf1[0][0]
    R1 = prf1[1][0]
    F11 = prf1[2][0]
    tp1 = CM1[0][0]
    fn1 = CM1[0][1]
    fp1 = CM1[1][0]
    tn1 = CM1[1][1]
    measures = [fp, assay, roc_auc1, MCC, kappa1, P1, R1, F11, A_count1, N_count1, tp1, fn1, fp1, tn1]
    with open(outdir, 'a') as outf:
        outf.write('{}'.format(m))
        for val in measures:
            outf.write('\t{}'.format(val))
        outf.write('\n')

#import os
import time
import pandas as pd
from sklearn import metrics

start = time.time()

#input data dir
file_path = 'C:/CESFP_project/CrossValidation/Assay_{}/{}_{}_full_predictions.txt'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
#sel_assays = ['527']
outdir = 'CrossValidation/CM.csv'

##load list of assays
#assay_list=[]
#with open(sel8_IDs, 'r') as f:
#    for line in f:
#        assay_list.append(line.strip())
#print(assay_list)

fp_types = ['htsfp', 'ecfp', 'cesfp']
#colors_list = ['orange','royalblue','forestgreen']
m = 0    
#Metric_dict = {}
headers = ['idx','fptype', 'assay', 'roc_auc', 'MCC', 'kappa', 'Precision', 'Recall', 'F1', 'A_count', 'N_count', 'tp', 'fn', 'fp', 'tn']
with open(outdir, 'w') as outf:
    outf.write('{}\n'.format('\t'.join(headers)))
for fp in fp_types:    
    print('Doing FP: {}'.format(fp))
    for assay in sel_assays: 
        print('Assay: {}'.format(assay))
        do_maths(assay, fp, m)
        
#        Metric_dict[m] = [fp, assay, j, roc_auc1, MCC, kappa1, P1, R1, F11, A_count1, N_count1, tp1, fn1, fp1, tn1]
        m += 1
    
#print('saving dataframes containing all scores')
#df = pd.DataFrame.from_dict(Metric_dict, orient='index')
#df.columns = ['fptype', 'assay', 'CV run','roc_auc', 'MCC', 'kappa', 'Precision', 'Recall', 'F1', 'A_count', 'N_count', 'tp', 'fn', 'fp', 'tn']
##dfav = df.mean()
#df.to_csv(outdir+'All_metrics.csv', sep='\t')
#dfav.to_csv(outdir+'{}_{}_metrics_average.csv'.format(fp,assay,j), sep='\t')



