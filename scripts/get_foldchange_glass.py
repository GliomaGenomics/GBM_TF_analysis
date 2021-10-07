# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:50:55 2020

@author: medgnt
"""

import pandas as pd
import numpy as np

data = pd.read_table('glass_data/gene_tpm_matrix_all_samples.tsv',sep='\t',header=0,index_col=0)
patients= pd.read_table('patient_lists/glass_gbm_idhwt_or_unknown_notstead_rna.txt',sep='\t', header=None)
ens=pd.read_table('downloaded_data/ensembl_v75_geneidtoname.txt',sep='\t', header=None)

#read in dictionary of gene to ensids
ensd={}
for r, en in ens.iterrows():
    ensd[en[1]]=en[0][0:15]

#convert data indexes from genenames to ensids
data=data.rename(ensd,axis='rows')

#keep only patients of interest
for c in data.columns:
    de=True
    for p in patients[0]:
        p=p.replace('-','.')
        if p+'.TP' in c or p+'.R1' in c:
            de=False
    if de==True:
        del data[c]

#filter genes
all_values=[]
for col in data:
    all_values.extend(list(data[col]))
all_values=pd.DataFrame(all_values)
lq=all_values[all_values >0].quantile(0.25)
print("lq="+str(lq))
ps=[]
for i in range(0,len(data)):
    abovep=data.iloc[i][data.iloc[i] >= lq[0]][data.iloc[i][data.iloc[i] >= lq[0]].index.str.contains('.TP')].count()
    totalp=data.iloc[i][data.iloc[i].index.str.contains('.TP')].count()
    proportionp=abovep/totalp
    abover=data.iloc[i][data.iloc[i] >= lq[0]][data.iloc[i][data.iloc[i] >= lq[0]].index.str.contains('.R1')].count()
    totalr=data.iloc[i][data.iloc[i].index.str.contains('.R1')].count()
    proportionr=abover/totalr
    if proportionp >=0.2 or proportionr >=0.2:
        ps.append('True_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
    else:
        ps.append('False_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))

data['ProportionAboveLQ']=ps
data=data[data['ProportionAboveLQ'].str.contains('True')]


for patient in patients[0]:
    patient=patient.replace('-','.')
    prim=data[data.columns[data.columns.str.startswith(patient+'.TP')][0]]
    recu=data[data.columns[data.columns.str.startswith(patient+'.R1')][0]]
    absolute=pd.DataFrame()
    actual=pd.DataFrame()
    absolute['genes']=data.index.str.split('.').str[0]
    actual['genes']=data.index.str.split('.').str[0]
    absolute['values']=list(abs(np.log2((recu+0.01)/(prim+0.01))))
    actual['values']=list(np.log2((recu+0.01)/(prim+0.01)))
    absolute.to_csv('ranks/glass_absolute_log2fc/'+patient.replace('.','-')+'.rnk',sep='\t',header=False,index=False)
    actual.to_csv('ranks/glass_actual_log2fc/'+patient+'.rnk',sep='\t',header=False,index=False)

