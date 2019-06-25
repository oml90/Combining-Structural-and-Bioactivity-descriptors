# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 16:55:55 2018

@author: kncv078

Input refined data from Noé script and annotation data (cmpd counts per assay). to make htsfp. also define threshold for minimum nuber of compounds per assay. 
"""

import pandas as pd
import time
import argparse
#import csv
import os
#from scipy import sparse
#import pubchempy as pcp

print('---'*30)
start = time.time()

# In[1] define inputs and output paths
ap = argparse.ArgumentParser(description='Takes refined hts data and generates a flag based fingerprint')
ap.add_argument('-a', action="store", dest='annotation_data', required=True, help='Path to annotation file')    # file with counts of compounds per assay post refinement with Noés script
ap.add_argument('-i', action="store", dest='hts_data', required=True, help='Path to refined hts data')          # folder with refined data produced by Noés script.
ap.add_argument('-o', action="store", dest='outpath', default='htsfp_gen_out/' , help='Optional: set output path for fingerprint file')
ap.add_argument('-t', action="store", dest='threshold', default=20000, type=int, help='Optional: threshold for minimum number of compounds per assay')

args = ap.parse_args()

if args.hts_data[-1] != '/':
    args.hts_data = args.hts_data+'/'

# additional inputs: file created by script: pubchem website or pubchempy
cid2smi_file    = 'C:/CESFP_project/CID2smi_pubchem.txt'

# set path for output
#if not os.path.exists(args.outpath):
#    os.makedirs(args.outpath)
outpath = args.outpath+'htsfp_t'+str(args.threshold)

# In[2] Load files
# load Annotation file and select assays with more than N cmpds using info from annotation file
print('loading annotation file...')
statsdf = pd.read_csv(args.annotation_data, sep='\t')
#statsdf.columns = ['assay','count']
#assay_selection = statsdf.loc[statsdf['count']>=args.threshold]['assay'].values
assay_selection = list(statsdf['assay'][statsdf['count']>=args.threshold])
print('Number of assays containing at least', args.threshold, 'Assays:', len(assay_selection))
print(time.time() - start)

# load compoundID 2 canonical smiles file 
print('loading CID 2 smiles file...')
cid2smi = {}
with open(cid2smi_file, 'r') as f:
    for line in f:
        linels = line.strip().split('\t')
        cid2smi[linels[0]] = linels[1]
print(time.time() - start)

# In[4] MAIN PART: Make dictionary containing cmpd_ID : list(assay activity flags)
#   PART 1/2: make empty dict with all values x's of correct dimensions
print('Creating empty dict...')
main_cmpd_dict = {}
for cmpd in cid2smi:
    main_cmpd_dict[cmpd] = ['x']*len(assay_selection)
#    main_cmpd_dict[cmpd] = ['0']*len(assay_selection)
print('done.')    
print('dimensions of htsfp:\t{} x {}'.format(len(main_cmpd_dict),len(assay_selection)))
print(time.time() - start)

# In[5] PART 2/2: fill dictionary with 'flags' from refined HTS data
print('Filling dict with flag values...')
idx = -1
Alength = len(assay_selection)
miss_cid = 0
for assay in assay_selection:
    idx += 1
    if idx%10==0:print('... {}%'.format(round((idx+1)/Alength*100),0), end='\r')
    assay_file = args.hts_data+assay
    with open(assay_file, 'r') as f:
        j = 0
        for line in f:
            j += 1
            if j == 1: continue #header row
            line_list = line.strip().split('\t')
            Cid = line_list[0]
            flag = line_list[1]
            if flag not in ['A','N']: print('{} unknown flag: {}'.format(assay, flag), end='\n'); continue
            if Cid == '': miss_cid += 1; continue
            elif Cid in main_cmpd_dict:
                main_cmpd_dict[Cid][idx] = flag
#                if flag == 'A': main_cmpd_dict[Cid][idx] = '1'
            else:
                print('{}, CID not in main cmpd dict: {}'.format(assay, Cid))
                continue
        
print('number of entries without CID: {}'.format(miss_cid))
print('\ndone')
print(time.time() - start)

# In[6] save output file containing: cmpd ID: HTSFP
print('saving output files')

with open(outpath+'.txt', 'w') as csv_file:
    for key, value in sorted(main_cmpd_dict.items()):
       csv_file.write('{} {}\n'.format(key,' '.join(value)))
print('HTSFP matrix file saved to:', outpath+'.txt')
       
# save a file containing an ordered list of assay IDs used to generate this HTSFP
print('saving assay list to file')
with open(outpath+'-Assay_list.txt','w') as assay_list:
    for item in assay_selection:
        assay_id = item.split('.')[0]
        assay_list.write("{}\n".format(assay_id))
print('Assay list saved to:', outpath+'-Assay_list.txt')

print('Total time taken:',int(time.time() - start),'seconds')    
