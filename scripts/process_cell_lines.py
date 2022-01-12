import pandas as pd
import numpy as np

mrna_gene_list=pd.read_table('ranks/filtered_genelist_mrna.txt',sep='\t',header=None, index_col=0)

if False:
    for name in ['A172','GBM63']:
        data = pd.read_table('cell_lines/raw_data/'+name+'_CON_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
        data.index=[i[:15] for i in data.index]
        output=pd.DataFrame(index=data.index)
        for ii in [0,1,2]:
            prim=data[name+'_U_'+str(ii)]
            recu=data[name+'_T_'+str(ii)]
            output[name+'_'+str(ii)]=list(np.log2((recu+0.01)/(prim+0.01)))
            output[name+'_'+str(ii)].to_csv('cell_lines/ranks/actual_log2fc/'+name+'_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
            abs(output[name+'_'+str(ii)]).to_csv('cell_lines/ranks/absolute_log2fc/'+name+'_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
        output.index=data.index
        output=pd.merge(output,mrna_gene_list, how = "inner",left_index=True, right_index=True)
        output.to_csv('cell_lines/tables/log2fc_'+name+'.txt',sep='\t',header=True,index=True)

if True:    
    data = pd.read_table('cell_lines/raw_data/GBM63_CUTRUN2021_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    output=pd.DataFrame(index=data.index)
    for ii in [0,1]:
        prim=data['U_'+str(ii)]
        recu=data['T_'+str(ii)]
        output['GBM63_CUTRUN_'+str(ii)]=list(np.log2((recu+0.01)/(prim+0.01)))
        output['GBM63_CUTRUN_'+str(ii)].to_csv('cell_lines/ranks/actual_log2fc/GBM63_CUTRUN_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
        abs(output['GBM63_CUTRUN_'+str(ii)]).to_csv('cell_lines/ranks/absolute_log2fc/GBM63_CUTRUN_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
    output.index=data.index
    output=pd.merge(output,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    output.to_csv('cell_lines/tables/log2fc_GBM63_CUTRUN.txt',sep='\t',header=True,index=True)
