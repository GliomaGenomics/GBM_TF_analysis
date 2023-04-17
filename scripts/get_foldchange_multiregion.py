# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 11:50:55 2020

@author: medgnt
"""

import pandas as pd
import numpy as np

data = pd.read_table('original_data/MultiRegion_all_geneFPKM.txt',sep='\t',header=0,index_col=0)
total_gene_list=pd.read_table('ranks/filtered_genelist_total.txt',sep='\t',header=None, index_col=0)
pairs = pd.read_table('original_data/MultiRegion_pairs.txt',sep='\t',header=0,index_col=0)

data.index=data.index.str.slice(stop=15)
data=pd.merge(data,total_gene_list, how = "inner",left_index=True, right_index=True)

data_out = pd.DataFrame(index=list(data.index))

for p in pairs.index:
    print(p)
    data_out[p]=list(np.log2((data[pairs.loc[p,'recurrent']]+0.01)/(data[pairs.loc[p,'primary']]+0.01)))
    abs(data_out[p]).to_csv('ranks/absolute_log2fc/'+p+'.rnk',sep='\t',header=False,index=True)
    data_out[p].to_csv('ranks/actual_log2fc/'+p+'.rnk',sep='\t',header=False,index=True)
 
