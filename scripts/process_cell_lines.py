import pandas as pd
import numpy as np

mrna_gene_list=pd.read_table('ranks/filtered_genelist_mrna.txt',sep='\t',header=None, index_col=0)

if False:
    for name in ['A172','GBM63']:
        data = pd.read_table('cell_lines/raw_data/'+name+'_CON_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
        data.index=[i[:15] for i in data.index]
        log2fc=pd.DataFrame(index=data.index)
        primary=pd.DataFrame(index=data.index)
        recurrent=pd.DataFrame(index=data.index)
        for ii in [0,1,2]:
            prim=data[name+'_U_'+str(ii)]
            recu=data[name+'_T_'+str(ii)]
            log2fc[name+'_'+str(ii)]=list(np.log2((recu+0.01)/(prim+0.01)))
            primary[name+'_'+str(ii)]=list(prim)
            recurrent[name+'_'+str(ii)]=list(recu)
            log2fc[name+'_'+str(ii)].to_csv('cell_lines/ranks/actual_log2fc/'+name+'_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
            abs(log2fc[name+'_'+str(ii)]).to_csv('cell_lines/ranks/absolute_log2fc/'+name+'_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
        log2fc.index=data.index
        log2fc=pd.merge(log2fc,mrna_gene_list, how = "inner",left_index=True, right_index=True)
        log2fc.to_csv('cell_lines/tables/log2fc_'+name+'.txt',sep='\t',header=True,index=True)
        primary.index=data.index
        primary=pd.merge(primary,mrna_gene_list, how = "inner",left_index=True, right_index=True)
        primary.to_csv('cell_lines/tables/primary_'+name+'.txt',sep='\t',header=True,index=True)
        recurrent.index=data.index
        recurrent=pd.merge(recurrent,mrna_gene_list, how = "inner",left_index=True, right_index=True)
        recurrent.to_csv('cell_lines/tables/recurrent_'+name+'.txt',sep='\t',header=True,index=True)

if False:    
    data = pd.read_table('cell_lines/raw_data/GBM63_CUTRUN2021_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    log2fc=pd.DataFrame(index=data.index)
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1]:
        prim=data['U_'+str(ii)]
        recu=data['T_'+str(ii)]
        log2fc['GBM63_CUTRUN_'+str(ii)]=list(np.log2((recu+0.01)/(prim+0.01)))
        primary['GBM63_CUTRUN_'+str(ii)]=list(prim)
        recurrent['GBM63_CUTRUN_'+str(ii)]=list(recu)
        log2fc['GBM63_CUTRUN_'+str(ii)].to_csv('cell_lines/ranks/actual_log2fc/GBM63_CUTRUN_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
        abs(log2fc['GBM63_CUTRUN_'+str(ii)]).to_csv('cell_lines/ranks/absolute_log2fc/GBM63_CUTRUN_'+str(ii)+'.rnk',sep='\t',header=False,index=True)
    log2fc.index=data.index
    log2fc=pd.merge(log2fc,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    log2fc.to_csv('cell_lines/tables/log2fc_GBM63_CUTRUN.txt',sep='\t',header=True,index=True)
    primary.index=data.index
    primary=pd.merge(primary,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    primary.to_csv('cell_lines/tables/primary_GBM63_CUTRUN.txt',sep='\t',header=True,index=True)
    recurrent.index=data.index
    recurrent=pd.merge(recurrent,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    recurrent.to_csv('cell_lines/tables/recurrent_GBM63_CUTRUN.txt',sep='\t',header=True,index=True)



if False:    
    data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week1_FPKM.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    log2fc=pd.DataFrame(index=data.index)
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for i in ['1','3']:
        data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week'+i+'_FPKM.txt',sep='\t',header=0,index_col=0)
        prim=data.loc[:, data.columns.str.startswith('U')].mean(axis=1)
        recu=data.loc[:, data.columns.str.startswith('T')].mean(axis=1)
        log2fc['PDspheroids_'+i]=list(np.log2((recu+0.01)/(prim+0.01)))
        primary['PDspheroids_'+i]=list(prim)
        recurrent['PDspheroids_'+i]=list(recu)
        log2fc['PDspheroids_'+i].to_csv('cell_lines/ranks/actual_log2fc/PDspheroids_'+i+'.rnk',sep='\t',header=False,index=True)
        abs(log2fc['PDspheroids_'+i]).to_csv('cell_lines/ranks/absolute_log2fc/PDspheroids_'+i+'.rnk',sep='\t',header=False,index=True)
    log2fc=pd.merge(log2fc,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    log2fc.to_csv('cell_lines/tables/log2fc_PDspheroids.txt',sep='\t',header=True,index=True)
    primary=pd.merge(primary,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    primary.to_csv('cell_lines/tables/primary_PDspheroids.txt',sep='\t',header=True,index=True)
    recurrent=pd.merge(recurrent,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    recurrent.to_csv('cell_lines/tables/recurrent_PDspheroids.txt',sep='\t',header=True,index=True)

if True:
    data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    log2fc=pd.DataFrame(index=data.index)
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for i in ['48','72']:
        data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
        prim=data.loc[:, data.columns.str.startswith('ME'+i+'hrU')].mean(axis=1)
        recu=data.loc[:, data.columns.str.startswith('ME'+i+'hrT')].mean(axis=1)
        log2fc['ME_'+i]=list(np.log2((recu+0.01)/(prim+0.01)))
        primary['ME_'+i]=list(prim)
        recurrent['ME_'+i]=list(recu)
        log2fc['ME_'+i].to_csv('cell_lines/ranks/actual_log2fc/ME_'+i+'.rnk',sep='\t',header=False,index=True)
        abs(log2fc['ME_'+i]).to_csv('cell_lines/ranks/absolute_log2fc/ME_'+i+'.rnk',sep='\t',header=False,index=True)
    log2fc=pd.merge(log2fc,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    log2fc.to_csv('cell_lines/tables/log2fc_ME.txt',sep='\t',header=True,index=True)
    primary=pd.merge(primary,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    primary.to_csv('cell_lines/tables/primary_ME.txt',sep='\t',header=True,index=True)
    recurrent=pd.merge(recurrent,mrna_gene_list, how = "inner",left_index=True, right_index=True)
    recurrent.to_csv('cell_lines/tables/recurrent_ME.txt',sep='\t',header=True,index=True)
