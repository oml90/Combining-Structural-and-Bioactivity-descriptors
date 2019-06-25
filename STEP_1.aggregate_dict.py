#!/usr/bin/env python3

# Author : Noe Sturm
# Edited : Oliver Laufkötter 
# Date   : 2018-03-23
# Aim    : Aggregates duplicated compounds (tested at same concentration) of HTS assays and returns activity flags from Z-score or activity readout.
#
#	   1/ Determines the most common concentration of an HTS assay
#
#	   2/ Selects records with  most common concentration +/- tolerance (default tolerance=0.05)
#
#	   3/ Sets the "activity flag".
#			> Z_SCORE_RESULT_FLAG if defined
#			> RESULT_FLAG if Z_SCORE_RESULT_FLAG not defined.
#			> Discard record if none of the flags are defined
#
#	   4/ Aggregates flags of duplicated records: most common flag (A for active, N for negative).
#             		> A if count(N) = count(A).  
#
#          5/ Prints "pre-processed" HTS assay files for HTS-FP generation: 
#			> compound-id 
#			> activity flag
#			> concentration
#			> flag source
# 	   
#	  6/ Prints counts relative to assays in stdout:
#			> assay file 
#			> most common concentration
#			> tolerance (+/-)
#			> total number of record in the assay file
#			> number of records with most common concentration +/- tolerance
#			> number of unique compounds with most common concentration +/- tolerance
#			> time spent for processing the file


import time
import sys, os
#import numpy as np
#import pandas as pd

if len(sys.argv) != 3:
	print ('Usage : python ' + sys.argv[0] + ' [ Raw hts data  ] [ output folder ]')
	quit()

output = sys.argv[2]
if not os.path.isdir(output):
	print ('Usage : python ' + sys.argv[0] + ' [ Raw hts data  ] [ output folder ]')
	sys.stderr.write("ERROR: can't find output folder at {}\n".format(output))
	sys.exit()
if output[-1] != '/': output=output+'/'

# headers
header = "CID\tflag\n"


# READ the list of assay files
assay_list = []
for file in os.listdir(sys.argv[1]):
    assay_list.append(file)
#content = open(sys.argv[1]).read().splitlines()

start = time.time()
sys.stderr.write('assays:{}\n'.format(len(assay_list)))



# loop through the assays
K=0
for assay in assay_list:
    assay_path = sys.argv[1]+assay
    
    f_type = assay.split('.')[-1]
    f_string = list(assay)
    if f_type != 'tsv':
        print('file: {} ==> is not of type .tsv, moving to next file'.format(assay))
        continue
    K+=1
    sys.stderr.write('{}\t({}/{})\t\r'.format(assay,K,len(assay_list)))
    # check if file is already present: if present, can we skip => NO because the assay can be updated (there are some with more than one data).
    
#    if os.path.exists(assay_path) == True:
#            #print('Assay with name: {} already exists under this outpath, moving on to next assay. Change outpath or move/rename existing files'.format(assay))
#            continue
    if not os.path.isfile(assay_path):
        sys.stderr.write("ERROR: can't find file {}\n".format(assay_path))
        continue

    cmpd_data={}
    start_assay = time.time()
    assay_content = open(assay_path, encoding='latin-1').read().splitlines()
	
    n_cmpd_total = len(assay_content)
    n_cmpd_interest = 0
    first_line = False
    n_head_rows = 0
    for line in assay_content:
        if line[0] == '1': first_line = True; 
        if line[0] + line[1] == '1\t': print('\nnumber of head rows: {}'.format(n_head_rows))
        if first_line == True:
            linelist = line.split("\t")
#            Sid     = linelist[1]
            Cid     = linelist[2]
            flag    = linelist[3]
    		
            # no CID number => SKIP #use Sample ID: sid?
            if Cid=='': continue
            
            # FLAG SELECTION
            if flag == 'Inactive': flag = 'N'
            elif flag == 'Active': flag = 'A'
            elif flag == 'Inconclusive': continue
            elif flag == 'Unspecified': continue
            elif flag == '': print('{} No flag set for CID: {}'.format(assay, Cid)); continue
            else: print('{} unknown flag: {}'.format(assay, flag), end='\n'); continue
    
            n_cmpd_interest+=1		
 
            # APPEND COMPOUND RECORD(S) IN A DICTIONNARY
            if Cid not in cmpd_data: cmpd_data[Cid]={'flag':[],'assay':assay}
            cmpd_data[Cid]['flag'].append(flag)
        else: n_head_rows += 1
        
    # SKIP IF NO COMPOUNDS - should never happen
    if len(cmpd_data)==0:
        end_assay = time.time()
        sys.stdout.write(assay+"\t"+str(n_cmpd_total)+"\t"+str(n_cmpd_interest)+"\t"+str(len(cmpd_data))+"\t"+str(round(end_assay-start_assay,5))+"\n")
        continue

    # AGGREGATION
    # flag : most common value, if equal number of As and Ns, select A.
    fout=open(output+assay.split('/')[-1],'w')
    fout.write(header)
    for cmpd in cmpd_data:
        # compound has multiple records
        if len(cmpd_data[cmpd]['flag']) > 1:
            flag = '//'
            n = [x for x in cmpd_data[cmpd]['flag'] if x=='N']
            if len(n) > len(cmpd_data[cmpd]['flag'])/2: flag = 'N' # if there are more Ns than half of the total number of records for the compound, select N.
            else:flag='A'
            fout.write(cmpd+"\t"+flag+"\n")
        # compound has single record
        else:	
            fout.write(cmpd+"\t"+cmpd_data[cmpd]['flag'][0]+"\n")
    fout.close()
	
    # log
    end_assay = time.time()
    sys.stdout.write(assay+"\t"+str(n_cmpd_total)+"\t"+str(n_cmpd_interest)+"\t"+str(len(cmpd_data))+"\t"+str(round(end_assay-start_assay,5))+"\n")

# save compound set to file
print("Elapsed time: {}".format(time.time()-start))


