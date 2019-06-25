# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 14:19:42 2018

@author: kncv078
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.sparse import coo_matrix, save_npz
import time
import numpy as np
#import my_functions.py

# In[1) FUNCTION: convert full htsfp cmpd list to rdkit ecfp, inputs are 1: smiles list, 2: ecfp radius  
def ecfp_from_smiles(cid_list, id2smiles_dict, ecfpRadius):
    print('\ngenerating ecfp from smiles...')
         
    smi2ecfp = {}
    bad_smiles = {}
    i = 0
    l = len(id2smiles_dict)
    for cmpd in id2smiles_dict:
        mol = Chem.MolFromSmiles(id2smiles_dict[cmpd])
        i += 1
        if not mol: #shouldnt be any, already filtered out in other step
            print(i, cmpd)
            bad_smiles[cmpd]=id2smiles_dict[cmpd]
        else:
            ecfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, radius=ecfpRadius, nBits=ecfplen))
            smi2ecfp[cmpd] = ecfp
        print('working: {}%'.format(round(i/l*100,1)),end='\r')
    print('\nNumber of bad smiles found:\t{}'.format(len(bad_smiles)))
    return (smi2ecfp, bad_smiles)


# In[2] Start, define starting time
start = time.time()
print('---'*30)

# In[3] INPUT VARIABLES: define input file path and ecfp radius (usually 2 or 3)
htsfp_cmpds_file = 'C:/CESFP_project/htsfp_sparse_cmpdlist.txt'
cid2smi_file = 'C:/CESFP_project/CID2smi_pubchem.txt'

ecfpRadius = 2 # 2=ecfp4 3=ecfp6
ecfplen = 1024
   
# In[4] load in htsfp column 1 for list of valid smiles
print('loading htsfp CID list...')
htsfp_CIDs = []
with open(htsfp_cmpds_file, 'r') as f:
    for line in f:
        htsfp_CIDs.append(line.strip())
print('Number of AZIDs:\t{}'.format(len(htsfp_CIDs)))

#Load AZ to canonical smiles dict
cid2smi = {}
print('loading cid 2 smiles dict')
with open(cid2smi_file, 'r') as f:
    for line in f:
        linels = line.strip().split('\t')
        cid2smi[linels[0]] = linels[1]
print('length of cid2smi file:\t{}'.format(len(cid2smi)))

# In[4] make sparse matrix of ECFPs
print('\ngenerating ecfp from smiles...')
row=[]
col=[]
data=[]
cmpdlist=[]

bad_smiles = {}
l = len(cid2smi)
for i,cmpd in enumerate(sorted(cid2smi)):
    mol = Chem.MolFromSmiles(cid2smi[cmpd])
    if not mol: #shouldnt be any, already filtered out in other step
        print('WARNING: bad smiles:',i, cmpd, cid2smi[cmpd])
        bad_smiles[cmpd]=cid2smi[cmpd]
    else:
        ecfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, radius=ecfpRadius, nBits=ecfplen))
    if i%20000==0: print('{}%'.format(round(i/l*100,0)),end='\r')
    cmpdlist.append(cmpd)
    for j, bit in enumerate(ecfp):
        if bit==1:
            row.append(i)
            col.append(j)
            data.append(np.int8(bit))

print('\nConverting to sparse matrix (coo)...')
cooM = coo_matrix((data, (row, col)), shape=(l, ecfplen))
save_npz('ecfp_sparse.npz',cooM)

#save corresponding cmpd list
with open('ecfp_sparse_cmpdlist.txt','w') as outf:
    for cmpd in cmpdlist:
        outf.write('{}\n'.format(cmpd))     
    
print('\nNumber of bad smiles found:\t{}'.format(len(bad_smiles)))

# In[5] run function: make ecfp from smiles list, returns dict.
#print('\ngenerating ecfp from smiles...')
#     
#bad_smiles = {}
#i = 0
#l = len(cid2smi)
#with open('Step3_out-cid2ecfp{}.txt'.format(ecfpRadius*2), 'w') as f:
#    for cmpd in cid2smi:
#        mol = Chem.MolFromSmiles(cid2smi[cmpd])
#        i += 1
#        if not mol: #shouldnt be any, already filtered out in other step
#            print(i, cmpd)
#            bad_smiles[cmpd]=cid2smi[cmpd]
#        else:
#            ecfp = list(AllChem.GetMorganFingerprintAsBitVect(mol, radius=ecfpRadius, nBits=1024))
#            f.write('{} {}\n'.format(cmpd, ' '.join(str(b) for b in ecfp)))
#        if i%10000==0: print('working: {}%'.format(round(i/l*100,1)),end='\r')
#print('\nNumber of bad smiles found:\t{}'.format(len(bad_smiles)))

# In[6] save ECFPs to file (non-sparse format - dense) very inefficient. but im not sure how to use sparse... yet.
#print('saving output files: ecfp file and bad smiles (if any)')
#with open('Step3_out-cid2ecfp{}.txt'.format(ecfpRadius*2), 'w') as f:
#    for key, value in cid2ecfp.items():
#        val = ''.join([str(x) for x in value])
#        f.write('{},{}\n'.format(key, val))
    
# save bad smiles to file; if any       
if len(bad_smiles) > 0:
    with open('bad_smiles_from_make_ecfp.txt', 'w') as f:
        for cid, smi in bad_smiles.items():
            f.write('{}\t{}\n'.format(cid, smi))
                 
print('done\ntime taken:\t{}'.format(time.time()-start))

