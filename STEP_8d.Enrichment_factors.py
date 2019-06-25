# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 14:13:21 2018

@author: kncv078
"""

#import os
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib as mpl
#mpl.use('Agg') #
#import matplotlib.pyplot as plt
#from sklearn import metrics
#from matplotlib_venn import venn2, venn3
#import seaborn
#from rdkit import Chem
#from rdkit.Chem.Scaffolds import MurckoScaffold
#from scipy.sparse import coo_matrix

start = time.time()
FP_sel = str(0) # 0 or 75
EFx_list = [0.002,0.003,0.005,0.0075,0.01,0.02,0.03,0.05]# enrichment factor proportion
EFx_list = [0.01]
cvflds = 6
#input data dir
#pred_files = '//samba.scp.astrazeneca.net/cc/kncv078/MESFP_projects/pubchem_dump_analysis/CrossValidation/crossValidation_RF_{}_analysis_sel8_Assays/'
pred_files = 'C:/CESFP_project/CrossValidation/Assay_{}/{}_{}_cvfold{}_predictions.txt'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
#cid2smi_file = 'C:/CESFP_project/CID2smi_pubchem.txt'

##load list of assays
#print('loading assay list..')
#assay_list=[]
#with open(sel8Assays_IDs_file, 'r') as f:
#    for line in f:
#        assay_list.append(line.strip())

# In[3] load all RF predictions for each of 9 new assays

labels1 = ['TP', 'TN', 'FP', 'FN', 'Recall', 'Precision', 'Accuracy', 'HR_all']
labels2 = ['EFx %','lenTopX', 'TP_2', 'FP_2', 'HR_topX', 'HR_all', 'EF']
ecfp_results = pd.DataFrame({'metric' : labels2})
htsfp_results  = pd.DataFrame({'metric' : labels2})
mesfp_results  = pd.DataFrame({'metric' : labels2})

print('beginning loop')
EF_dict = {}
for EFx in EFx_list: 
    print('---'*15)
    print('analysis of EF = {}%'.format(EFx*100))
    for assay in sel_assays:
        print('Doing assay: {}'.format(assay))
        FP_list = ['ecfp', 'htsfp', 'cesfp']
        for fp_name in FP_list:
            for CV_run in range(1,cvflds+1):
    
                FP_df = pd.read_csv(pred_files.format(assay,assay,fp_name,CV_run), sep='\t').sort_values('PredProb(A)', ascending=False).reset_index(drop=True)
            
                TP = len(FP_df[FP_df['label']=='A']) 
                count_all = len(FP_df)
                HR_all = TP/count_all
                        
                lenTopX = int(count_all*EFx) 
                Topx_df = FP_df.iloc[:lenTopX,:]
                TP_2 = len(Topx_df[Topx_df['label']=='A'])
                HR_topX = TP_2/lenTopX
                EF = HR_topX/HR_all
                
                EF_dict['{}_{}-{}_{}'.format(fp_name, assay, CV_run, (EFx*100))] = [fp_name, assay, EF, HR_topX, TP_2, lenTopX, HR_all, TP, count_all]
        
    
#print('{} \n {} \n {}'.format(ecfp_df, htsfp_df, mesfp_df))
#ecfp_df = ecfp_df.set_index('metric')
#htsfp_df = htsfp_df.set_index('metric')
#mesfp_df = mesfp_df.set_index('metric')


#ecfp_df.to_csv('/Random_Forest_hts0/Enrichment_ecfp.txt', sep='\t')
#htsfp_df.to_csv('/Random_Forest_hts0/Enrichment_htsfp.txt', sep='\t')
#mesfp_df.to_csv('/Random_Forest_hts0/Enrichment_mesfp.txt', sep='\t')

EF_df = pd.DataFrame.from_dict(EF_dict, orient='index')
EF_df.columns = ['fptype','AID','EF', 'HR_topX', 'TP_topX', 'count_topX', 'HR_all', 'TP_all', 'count_all']

d = {}
for assay in sel_assays:
    df = EF_df[EF_df['AID']==assay]
    for fpt in ['htsfp','ecfp','cesfp']:
        vals = list(df['EF'][df['fptype']==fpt])
        mean = np.mean(vals)
        stddev = np.std(vals)
        d[assay+fpt] = [assay,fpt,mean,stddev]
df=pd.DataFrame.from_dict(d,orient='index')
df.columns = ['assay','fptype','EF_mean','EF_sd']


#EF_averages = pd.DataFrame(columns=['EF', 'HR_topX', 'TP_topX', 'count_topX', 'HR_all', 'TP_all', 'count_all'])
#for i in range(0,1920,6):
##    print(i)
#    averages = EF_df.iloc[i:i+6,:].mean().values
#    fp_type = FP_list[i%3]
#    assay = sel_assays[i//30%8]
#    EFx = EFx_list[i//240]
#    EF_averages.loc['{}_{}_{}_mean'.format(fp_type, assay, EFx*100)] = averages
    

#EF_df.to_csv('CrossValidation/Enrichment_factors.csv' ,sep='\t')
#EF_averages.to_csv('CrossValidation/Enrichment_factor_averages.csv' ,sep='\t')


ecMeans = list(df['EF_mean'][df['fptype']=='ecfp'])
ecStd = list(df['EF_sd'][df['fptype']=='ecfp'])
htsMeans = list(df['EF_mean'][df['fptype']=='htsfp'])
htsStd = list(df['EF_sd'][df['fptype']=='htsfp'])
bashMeans = list(df['EF_mean'][df['fptype']=='cesfp'])
bashStd = list(df['EF_sd'][df['fptype']=='cesfp'])

plt.figure(figsize=[18,6])
N = len(ecMeans)


ind = np.arange(N)    # the x locations for the groups
width = 0.25         # the width of the bars

p1 = plt.bar(ind, ecMeans, width, bottom=0, yerr=ecStd)
p2 = plt.bar(ind + width, htsMeans, width, bottom=0, yerr=htsStd)
p3 = plt.bar(ind + width*2, bashMeans, width, bottom=0, yerr=bashStd)

#plt.title('ROC-AUC', size=18)
plt.xticks(ind + width, sel_assays)

plt.legend((p1[0], p2[0], p3[0]), ('ECFP', 'HTSFP', 'BaSH'), prop={'size': 18})

plt.xlabel('assay ID', size=22)          #(HTSFP: 0-553, ECFP4: 554-1577)
plt.ylabel('Enrichment Factor {}%'.format(EFx_list[0]*100), size=22)
plt.ylim(0)
#plt.xlim([-6,tlen+6])
plt.tick_params(axis='both', labelsize=18)
plt.tick_params(axis='x', labelrotation=45)

plt.savefig('EF{}%_barplot.png'.format(EFx_list[0]*100), bbox_inches='tight')
