# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 20:15:35 2019

@author: laufkoetter
"""

# plot the performance metrics into beautiful bar pplots

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

results_file = 'C:/CESFP_project/CrossValidation/All_metrics.csv'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
#sel_assays = ['527', '798', '2606', '463104','504812', '588497','623901', '686964']

df = pd.read_csv(results_file, sep='\t').set_index('idx')


d = {}
for assay in sel_assays:
    df2 = df[df['assay']==int(assay)]
    for fpt in ['htsfp','ecfp','cesfp']:
        d[assay+fpt] = [assay,fpt]
        for m in ['roc_auc','MCC','kappa','Precision','Recall','F1']:
            vals = list(df2[m][df2['fptype']==fpt])
            mean = np.mean(vals)
            stddev = np.std(vals)
            d[assay+fpt].append(mean)
            d[assay+fpt].append(stddev)
df=pd.DataFrame.from_dict(d,orient='index')
df.columns = ['assay','fptype','ra_mean','ra_sd','MCC_mean','MCC_sd','k_mean','k_sd','P_mean','P_sd','R_mean','R_sd','F1_mean','F1_sd']

measures = {'ra_':'ROC-AUC','MCC_':'MCC','k_':'Kappa-score','P_':'Precision','R_':'Recall','F1_':'F1-score'}
for met in measures:
    ecMeans = list(df[met+'mean'][df['fptype']=='ecfp'])
    ecStd = list(df[met+'sd'][df['fptype']=='ecfp'])
    htsMeans = list(df[met+'mean'][df['fptype']=='htsfp'])
    htsStd = list(df[met+'sd'][df['fptype']=='htsfp'])
    bashMeans = list(df[met+'mean'][df['fptype']=='cesfp'])
    bashStd = list(df[met+'sd'][df['fptype']=='cesfp'])
    
    plt.figure(figsize=[18,6])
    N = len(ecMeans)
    
    
    ind = np.arange(N)    # the x locations for the groups
    width = 0.25         # the width of the bars
    
    p1 = plt.bar(ind, ecMeans, width, bottom=0, yerr=ecStd)
    p2 = plt.bar(ind + width, htsMeans, width, bottom=0, yerr=htsStd)
    p3 = plt.bar(ind + width*2, bashMeans, width, bottom=0, yerr=bashStd)
    
    #plt.title('ROC-AUC', size=18)
    plt.xticks(ind + width, sel_assays)
    #plt.xticklabels(sel_assays)
    
    plt.legend((p1[0], p2[0], p3[0]), ('ECFP', 'HTSFP', 'BaSH'), prop={'size': 18})
    
    plt.xlabel('assay ID', size=22)          #(HTSFP: 0-553, ECFP4: 554-1577)
    plt.ylabel(measures[met], size=22)
    #plt.xticks(range(0,1600,100))
    #plt.plot(range(htslen,tlen),fimp_max[htslen:], label='ECFP4 max', color='darkblue')
    #plt.plot(range(htslen,tlen),fimp_av[htslen:], label='ECFP4 mean', color='deepskyblue')
    #plt.plot(fimp_max[:htslen], label='HTSFP max', color='orangered')
    #plt.plot(fimp_av[:htslen], label='HTSFP mean', color='sandybrown')
    #plt.legend(loc="upper right", ncol=2, prop={'size': 14}, frameon=True, facecolor='gainsboro') 
    if met == 'ra_':
        plt.ylim(0.5,1.02)
    else:
        plt.ylim(0)
    #plt.xlim([-6,tlen+6])
    plt.tick_params(axis='both', labelsize=18)
    plt.tick_params(axis='x', labelrotation=45)
    #plt.grid(axis='y', color='silver', alpha=0.5)
    #plt.tick_params(axis='both', which='major', length=5, width=1.7)
    plt.savefig('{}_barplot.png'.format(measures[met]), bbox_inches='tight')