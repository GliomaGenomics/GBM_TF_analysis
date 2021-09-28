# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:50:55 2020

@author: medgnt
"""

import pandas as pd
import numpy as np

data = pd.read_table('original_data/PvR_isoformfpkm_all_LS_23062021.txt.txt',sep='\t',header=0,index_col=0)
met= pd.read_table('original_data/MetaData_LS_230621.txt',sep='\t',header=0,index_col=2)
total_gene_list=pd.read_table('ranks/filtered_tsslist_total.txt',sep='\t',header=None, index_col=0)
mrna_gene_list=pd.read_table('ranks/filtered_tsslist_mrna.txt',sep='\t',header=None, index_col=0)

tssi={}
with open('intermediate_files/transcript_to_tss_position.txt','r+') as file:
    for line in file:
        l=line.strip().split()
        tssi[l[0]]=l[1]


tss=[]
dro=[]
for t in data.index:
    if t in tssi:
        tss.append(tssi[t])
#    elif data.loc[t,'GeneName'] in tssn:
#        tss.append(tssn[data.loc[t,'GeneName']])
    else:
        dro.append(t)
        print(t)

data.drop(dro,inplace=True)

data['tss']=tss
data2=data.groupby(['tss']).sum()
data=data2

data_total = data.copy(deep=True)
data_mrna = data.copy(deep=True)
for c in data.columns:
    id=c[:max([c.find('_R_'),c.find('_P_'),c.find('_Prim'),c.find('_Recu')])]
    if 'Total' in met['LibraryType'][id]:
        del data_mrna[c]
    elif 'mRNA' in met['LibraryType'][id]:
        del data_total[c]

d={}
d['total']=data_total
d['mrna']=data_mrna

for di in d:

    if di=='total':
        d[di]=pd.merge(d[di],total_gene_list, how = "inner",left_index=True, right_index=True)
    if di=='mrna':
        d[di]=pd.merge(d[di],mrna_gene_list, how = "inner",left_index=True, right_index=True)


    for patient in met.index:
        prim=''
        recu=''
        if patient +'_Primary_FPKM' in d[di].columns:
            prim=d[di][patient+'_Primary_FPKM']
            recu=d[di][patient+'_Recurrent_FPKM']
        elif patient +'_P_FPKM' in d[di].columns:
            prim=d[di][patient+'_P_FPKM']
            recu=d[di][patient+'_R_FPKM']
        else:
            continue
        absolute=pd.DataFrame()
        actual=pd.DataFrame()
        absolute['genes']=d[di].index
        actual['genes']=d[di].index
        absolute['values']=list(abs(np.log2((recu+0.01)/(prim+0.01))))
        actual['values']=list(np.log2((recu+0.01)/(prim+0.01)))
        absolute.to_csv('ranks/absolute_log2fc_tss_newrelease/'+patient+'.rnk',sep='\t',header=False,index=False)
        actual.to_csv('ranks/actual_log2fc_tss_newrelease/'+patient+'.rnk',sep='\t',header=False,index=False)
