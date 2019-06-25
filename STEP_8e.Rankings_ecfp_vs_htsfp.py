# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:02:22 2018

script to compare rankings of active compounds from the ECFP4 vs the HTSFP

@author: kncv078
"""

import pandas as pd
import matplotlib.pyplot as plt

top_Number = 1000

#input data dir
pred_files = 'C:/CESFP_project/CrossValidation/Assay_{}/{}_{}_full_predictions.txt'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']

#sel8_names = mf.file_lines_to_list(sel8_names_file)
#sel_assays=[int(x) for x in sel_assays]
#sel_assays.sort()

for assay in sel_assays:
    print('doing AID:{}'.format(assay))
    BaSH = pd.read_csv(pred_files.format(assay,assay,'cesfp'), sep='\t').sort_values('PredProb(A)', ascending=False).set_index('cmpd')
    ecfp = pd.read_csv(pred_files.format(assay,assay,'ecfp'), sep='\t').sort_values('PredProb(A)', ascending=False).set_index('cmpd')
    htsfp = pd.read_csv(pred_files.format(assay,assay,'htsfp'), sep='\t').sort_values('PredProb(A)', ascending=False).set_index('cmpd')
    
    BaSH['rank']=range(len(BaSH))
    ecfp['rank']=range(len(ecfp))
    htsfp['rank']=range(len(htsfp))
    
    sel_cmpds = BaSH[:top_Number]
    sel_A = sel_cmpds[sel_cmpds['label']=='A']
    sel_N = sel_cmpds[sel_cmpds['label']=='N']
     
    sC_ecfp_A = ecfp.loc[sel_A.index]
    sC_htsfp_A = htsfp.loc[sel_A.index]
    sC_ecfp_N = ecfp.loc[sel_N.index]
    sC_htsfp_N = htsfp.loc[sel_N.index]
#    lowest_rank = max([max(sC_ecfp_A['rank']),max(sC_ecfp_N['rank']),max(sC_htsfp_A['rank']),max(sC_htsfp_N['rank'])])
    lowest_rank = 500000
    
    plt.figure(figsize=(12,8))
    ax = plt.gca()
    plt.title('AID:{}'.format(assay), size=26)
    plt.plot([0,lowest_rank],[0,lowest_rank], color='k', linewidth=2, alpha=0.8)
    plt.plot([1000,lowest_rank+100000],[1000,1000], color='k', linewidth=3, linestyle='--', alpha=0.4)
    plt.plot([1000,1000],[1000,lowest_rank+100000], color='k', linewidth=3, linestyle='--', alpha=0.4)
    plt.scatter(sC_ecfp_N['rank'],sC_htsfp_N['rank'], label='Inactives', color='orange', alpha=0.6)
    plt.scatter(sC_ecfp_A['rank'],sC_htsfp_A['rank'], label='Actives', color='green', alpha=0.85)
    plt.xlim([0.9,lowest_rank+100000])
    plt.ylim([0.9,lowest_rank+100000]) 
    plt.legend(loc="lower right", ncol=1, prop={'size': 22}, frameon=True, facecolor='lightgrey')  
    plt.xlabel('ECFP4 rank', size=24)
    plt.ylabel('HTSFP rank', size=24)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.tick_params(axis='both', which='major', length=6, width=2.0)
    ax.tick_params(axis='both', which='minor', length=4, width=1.5)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.grid(axis='both')
    plt.savefig('Rankings_plot_{}.png'.format(assay), bbox_inches='tight')
    plt.show()
#    break