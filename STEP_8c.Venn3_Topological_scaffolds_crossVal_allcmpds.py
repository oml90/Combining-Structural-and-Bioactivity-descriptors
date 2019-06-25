# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 13:59:38 2018

@author: kncv078
"""

def load_results(pred_files, assay, prcnt):
    df_ecfp = pd.read_csv(pred_files.format(assay,assay,'ecfp'), sep='\t').sort_values('PredProb(A)', ascending=False).reset_index(drop=True)
    if prcnt <= 100:
        numba = int(len(df_ecfp)*(prcnt/100))
    else:   numba = prcnt
    TS_ec_df = df_ecfp.iloc[:numba,:].set_index('cmpd')
    df_htsfp = pd.read_csv(pred_files.format(assay,assay,'htsfp'), sep='\t').sort_values('PredProb(A)', ascending=False).reset_index(drop=True)
    TS_hts_df = df_htsfp.iloc[:numba,:].set_index('cmpd')
    df_cesfp = pd.read_csv(pred_files.format(assay,assay,'cesfp'), sep='\t').sort_values('PredProb(A)', ascending=False).reset_index(drop=True)
    TS_ces_df = df_cesfp.iloc[:numba,:].set_index('cmpd')
    return TS_ec_df, TS_hts_df, TS_ces_df, numba

    

import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold as ms

start = time.time()

topX = 1000 #top x percent or just top x
use_TPs = True

# In[1] define input data dir
pred_files = 'C:/CESFP_project/CrossValidation/Assay_{}/{}_{}_full_predictions.txt'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
cid2smi_file = 'C:/CESFP_project/CID2smi_pubchem.txt'

# In[2] load files
##load list of assays
#print('loading assay list..')
#assay_list=[]
#i = 0
#with open(sel8Assays_IDs_file, 'r') as f:
#    for line in f:
#        i += 1
#        assay_list.append(line.strip())

# load dict of Compound ID 2 smiles - required for RDkit to make scaffolds
print('loading smiles file...')
cid2smi = {}
i = 0
with open(cid2smi_file, 'r') as f:
    for line in f:
        i += 1
        cid, smi = line.strip().split('\t')
        cid2smi[cid] = smi

# In[3] load all RF predictions for each of 9 new assays
print('Top scoring {}% of predictions'.format(topX))
i = 1
d=[]
for assay in sel_assays:
    outpath = 'CrossValidation/Venn_diagrams_Topological_scafs_top_{}%_allcmpds/'.format(topX)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    #load files and select top x percent of hits to DF
    print('Doing assay: {}'.format(assay))
    print('\t# g.scafs   /   # CIDs'.format(topX))
    i += 1
    ecfp_r, htsfp_r, BaSH_r, numba = load_results(pred_files, assay, topX)
    
    #use only true positives
    if use_TPs == True:
        ecfp_r = ecfp_r[ecfp_r['label']=='A']
        htsfp_r = htsfp_r[htsfp_r['label']=='A']
        BaSH_r = BaSH_r[BaSH_r['label']=='A']
    
    # remove the CIDs of the top x compounds...
    TPec = list(ecfp_r.index)
    TPhts = list(htsfp_r.index)
    TPces = list(BaSH_r.index)
    
    # caluculate smiles and Topological scaffold for each compound
    # analysing topological scaffolds...
    cmpd_lists = {'ecfp':TPec, 'htsfp':TPhts, 'cesfp':TPces}
    Generic_sets_dict = {}
    for FP_name, cmpds in cmpd_lists.items():
        gen_scaf_set = set()
        for cid in cmpds:
            if str(cid) in cid2smi:
                m_scaf =  ms.GetScaffoldForMol(Chem.MolFromSmiles(cid2smi[str(cid)]))                    
                g_scaf = Chem.MolToSmiles(ms.MakeScaffoldGeneric(m_scaf))
                gen_scaf_set.add(g_scaf)
            else: 
                print('NA???, thats not meant to happen....')
                continue
        
        print('{}\t   {}     /    {}'.format(FP_name,len(gen_scaf_set),numba))
        Generic_sets_dict[FP_name] = gen_scaf_set
    
    ec_scafs = Generic_sets_dict['ecfp']
    hts_scafs = Generic_sets_dict['htsfp']
    ces_scafs = Generic_sets_dict['cesfp']
    
#    o = ec_scafs.union(hts_scafs)
#    l=len(o)
#    ol = len(o&ces_scafs)
#    cl = len(ces_scafs)
#    u = cl-ol
#    if cl == 0 or ol ==0: continue
#    d.append(u/l*100)

    plt.figure(i, figsize=(6,5))
    plt.title('AID:'+assay, size=18)

    venndiag = venn3([ec_scafs, hts_scafs, ces_scafs], set_labels = ('ECFP ({})'.format(len(ec_scafs)), 
                     'HTSFP ({})'.format(len(hts_scafs)), 'BaSH ({})'.format(len(ces_scafs))), alpha=0)
    vennC = venn3_circles([ec_scafs, hts_scafs, ces_scafs], linewidth=5)
    vennC[0].set_color('dodgerblue')
    vennC[1].set_color('orangered')
    vennC[2].set_color('limegreen')
    vennC[0].set_alpha(0.4)
    vennC[1].set_alpha(0.4)
    vennC[2].set_alpha(0.4)
    for text in venndiag.set_labels:
        if text == None: continue
        text.set_fontsize(15)
    for text in venndiag.subset_labels:
        if text == None: continue
        text.set_fontsize(15)

    
    plt.savefig('{}{}-{}%_VennDiagram.png'.format(outpath,assay,topX), bbox_inches='tight')
#    break
print('done')




