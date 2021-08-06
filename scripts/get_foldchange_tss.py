# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:50:55 2020

@author: medgnt
"""

import pandas as pd
import numpy as np

data = pd.read_table('original_data/PvR_isoformfpkm_all.txt',sep='\t',header=0,index_col=0)
met= pd.read_table('original_data/MetaData.txt',sep='\t',header=0,index_col=2)


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
    id=c[:c.find('_')]
    if 'Total' in met['LibraryType'][id]:
        del data_mrna[c]
    elif 'mRNA' in met['LibraryType'][id]:
        del data_total[c]

d={}
d['total']=data_total
d['mrna']=data_mrna

for di in d:
    all_values=[]
    for col in d[di]:
        all_values.extend(list(d[di][col]))
    all_values=pd.DataFrame(all_values)
    lq=all_values[all_values >1].quantile(0.25)
    ps=[]
    for i in range(0,len(d[di])):
        abovep=d[di].iloc[i][d[di].iloc[i] >= lq[0]][d[di].iloc[i][d[di].iloc[i] >= lq[0]].index.str.contains('_P')].count()
        totalp=d[di].iloc[i][d[di].iloc[i].index.str.contains('_P')].count()
        proportionp=abovep/totalp
        abover=d[di].iloc[i][d[di].iloc[i] >= lq[0]][d[di].iloc[i][d[di].iloc[i] >= lq[0]].index.str.contains('_R')].count()
        totalr=d[di].iloc[i][d[di].iloc[i].index.str.contains('_R')].count()
        proportionr=abover/totalr
        if proportionp >=0.5 or proportionr >=0.5:
            ps.append('True_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
        else:
            ps.append('False_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
    print(str(di)+"_"+str(lq[0]))
    data['ProportionAboveLQ_'+di]=ps

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

        data.index=absolute.index

        absolute=absolute[data['ProportionAboveLQ_'+di].str.contains('True')]
        actual=actual[data['ProportionAboveLQ_'+di].str.contains('True')]

        absolute.to_csv('ranks/absolute_log2fc_tss/'+patient+'.rnk',sep='\t',header=False,index=False)
        actual.to_csv('ranks/actual_log2fc_tss/'+patient+'.rnk',sep='\t',header=False,index=False)
#data.to_csv('PvR_genefpkm_all_filter_tss_stringent.txt',sep='\t',header=True,index=True)

