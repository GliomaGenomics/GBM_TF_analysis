import pandas as pd
import numpy as np

mrna_gene_list=pd.read_table('ranks/filtered_genelist_mrna.txt',sep='\t',header=None, index_col=0)

data = pd.read_table('cell_lines/raw_data/'+'A172'+'_CON_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
data.index=[i[:15] for i in data.index]
output=pd.DataFrame(index=data.index)
prim=data['A172_U_0']
recu=data['A172_U_1']
output['A172_control']=list(np.log2((recu+0.01)/(prim+0.01)))
output['A172_control'].to_csv('cell_lines/ranks/actual_log2fc/A172_control.rnk',sep='\t',header=False,index=True)
abs(output['A172_control']).to_csv('cell_lines/ranks/absolute_log2fc/A172_control.rnk',sep='\t',header=False,index=True)
    
