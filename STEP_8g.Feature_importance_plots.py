# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 11:51:34 2018
plot feature importance for the results of the random forest 
@author: kncv078
"""


#import my_functions_local as mf
#import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('default')

FI_files = 'C:/CESFP_project/CrossValidation/FI/Feature_importance_cesfp_{}_{}.out'
sel_assays = ['522', '527', '555', '560', '746', '798', '1006', '1273', '1515', '2129', '2280', '2540', '2544', '2553', 
              '2606', '463104', '504406', '504454', '504812', '588497', '602363', '623901', '624414', '686964', '720700']
HTSFP_names_file = 'C:/CESFP_project/htsfp_sparse_assaylist.txt'

plot = True
flds = 6
ecfplen = 1024

assay_list = []
with open(HTSFP_names_file, 'r') as inf:
    for line in inf:
        assay_list.append(line.strip())

htslen = len(assay_list)-1
tlen = htslen+ecfplen

for j, assay in enumerate(sel_assays):
    #    print('{}%'.format(j/len(test_labels)*100), end='\r')
    if plot == True:
        plt.figure(j,figsize=[10,5])
        plt.title('AID:{}'.format(assay), fontsize=20)
        plt.xlabel('FP bit number', size=18)          #(HTSFP: 0-553, ECFP4: 554-1577)
        plt.ylabel('Feature Importance', size=18)
        plt.xticks(range(0,1600,100))
#        plt.plot([553,553],[0,0.13], color='k', linestyle='--', alpha=0.4)
    
    for i in range(1,flds+1):
        FIs = [] #mf.file_lines_to_list(FI_files.format(assay,i))
        with open(FI_files.format(assay,i),'r') as inf:
            for line in inf:
                FIs.append(line.strip())
        if i == 1: fimp_max = FIs.copy() ; fimp_av = FIs.copy() ; #print('here1')
        else: 
            for k in range(len(FIs)):
                fimp_max[k] = max(float(FIs[k]),float(fimp_max[k]))
                fimp_av[k] = float(FIs[k])+float(fimp_av[k])
                if i == flds:
                    fimp_av[k] = fimp_av[k]/flds
#        plt.plot(range(553,1577),FIs[553:], label='ECFP4_max',alpha=0.3)
#        plt.plot(FIs[:553], label='HTSFP_max', alpha=0.3
    
    if plot == True:
        plt.plot(range(htslen,tlen),fimp_max[htslen:], label='ECFP4 max', color='darkblue')
        plt.plot(range(htslen,tlen),fimp_av[htslen:], label='ECFP4 mean', color='deepskyblue')
        plt.plot(fimp_max[:htslen], label='HTSFP max', color='orangered')
        plt.plot(fimp_av[:htslen], label='HTSFP mean', color='sandybrown')
        plt.legend(loc="upper right", ncol=2, prop={'size': 14}, frameon=True, facecolor='gainsboro') 
        plt.ylim([0,max(0.1,max(fimp_max)+0.01)])
        plt.xlim([-6,tlen+6])
        plt.tick_params(axis='both', labelsize=14)
        plt.tick_params(axis='x', labelrotation=45)
        plt.grid(axis='y', color='silver', alpha=0.5)
        plt.tick_params(axis='both', which='major', length=5, width=1.7)
#        plt.show()
        plt.savefig('AID{}_Feature_Importance.png'.format(assay), bbox_inches='tight')

    # print out information on top 5 most important assays.
    topFIi = sorted(range(len(fimp_av[:htslen])), key=lambda k: fimp_av[k], reverse=True)[:5]
    print('AID{}\tindex\tImp_value'.format(assay))
    adjassaylist = assay_list.copy() 
    adjassaylist.remove(assay)
    for h in topFIi:
        print('{}\t{}\t{}'.format(adjassaylist[h], h, round(fimp_av[h], 3)))
    with open('FI_top5.csv','a') as outf:
        outf.write('AID{}\tindex\tImp_value\n'.format(assay))
        for h in topFIi:
            outf.write('{}\t{}\t{}\n'.format(adjassaylist[h], h, round(fimp_av[h], 3)))
#    break
