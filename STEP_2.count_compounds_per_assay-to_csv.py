# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 10:31:48 2018

@author: kncv078
"""

# Script does a count of the unique compounds in each assay of a directory eg. full_dump or refined_dump
# 

import pandas as pd
import time
import os

# In[1]
# function for cid count in refined pubchem dump
def assay_uniq_cmpd_counts(folder):
    all_cid_uniq = set()
    cid_lc_dict = {}
    i = 0
    for file in os.listdir(folder):
        i += 1
        if file.split('.')[1] != 'tsv':
            print('bad file: {}'.format(file))
            continue
        with open(folder+file, 'r', encoding='latin_1') as f:
            cid_set = []
            for line in f:
                if line[0] == 'C': continue
                line_ls = line.split('\t')
                cid = line_ls[0]
                flag = line_ls[1]
                if cid == '' and flag != '': # shouldnt happen
                    print('ERROR: {} no cid for assay {}'.format(file, cid))
                    continue
                else: cid_set.append(cid); all_cid_uniq.add(cid)
                
        print('{} \tCid repeats: {}'.format(file, len(cid_set)-len(set(cid_set))))
        cid_lc_dict[file] = len(set(cid_set))
        #print('processing: {}%   AID: {}   cid count: {}   sid count: {}'.format(int(i/len(os.listdir(folder))*100), file, len(cid_set), len(sid_set)), end='\r')
    CID_uniq_df = pd.DataFrame.from_dict(cid_lc_dict, orient='index').sort_index()
    CID_uniq_df.columns = ['count']
    return CID_uniq_df, all_cid_uniq

# function for cid count in RAW pubchem dump
def assay_uniq_cmpd_counts_RAW(folder):
    all_cid_uniq = set()
    cid_lc_dict = {}
    i = 0
    for file in os.listdir(folder):
        i += 1
        if file.split('.')[1] != 'tsv':
            print('bad file: {}'.format(file))
            continue
        with open(folder+file, 'r', encoding='latin_1') as f:
            cid_set = []
            startline = False
            for line in f:
                if line[0] == '1': startline = True
                if startline == True:
                    cid = line.strip().split('\t')[2]
                    if cid == '': 
                        print('{} no cid for assay {}'.format(file, cid))
                        continue
                    else: all_cid_uniq.add(cid)
                
        print('{} \tCid repeats: {}'.format(file, len(cid_set)-len(set(cid_set))))
        cid_lc_dict[file] = len(all_cid_uniq)
        #print('processing: {}%   AID: {}   cid count: {}   sid count: {}'.format(int(i/len(os.listdir(folder))*100), file, len(cid_set), len(sid_set)), end='\r')
    CID_uniq_df = pd.DataFrame.from_dict(cid_lc_dict, orient='index').sort_index()
    CID_uniq_df.columns = ['count']
    return CID_uniq_df, all_cid_uniq


# In[3]
start = time.time()

pubchem_dump = 'C:/CESFP_project/Step1a_out_refined_HTS_data/'
#pubchem_dump = '//samba.scp.astrazeneca.net/kncv078/HTS_Project-Fingerprints/pubchemAssaysRAW/'

print('beginning analysis')
dfc, cid_list = assay_uniq_cmpd_counts(pubchem_dump)
#dfc, cid_list = assay_uniq_cmpd_counts_RAW(pubchem_dump)

print('d1 done, time taken: {}'.format(time.time()-start))

with open('Step1b_out-all_uniq_cmpds_list.txt','w') as f:
    for i in cid_list:
        f.write('{}\n'.format(i))       
        
dfc.to_csv('Step1b_out-annotationData_uniq_cid_count_perAssay.csv', sep='\t')

print('finished, time taken: {}'.format(time.time()-start))