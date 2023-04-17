import pandas as pd
import numpy as np

all_prim=pd.DataFrame()
all_recu=pd.DataFrame()

data = pd.read_table('cell_lines/raw_data/M059K_GBM63ser_UvTFPKM.txt',sep='\t',header=0,index_col=0)
data=data[data.index.str.contains('_PAR_Y')==False]
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for ii in [2,3]:
    primary['GBM63_ser_'+str(ii)]=data['GBM63_ser_Rep'+str(ii)+'_U']
    recurrent['GBM63_ser_'+str(ii)]=data['GBM63_ser_Rep'+str(ii)+'_T']
for ii in [1,2]:
    primary['M059K_'+str(ii)]=data['M059K_Rep'+str(ii)+'_U']
    recurrent['M059K'+'_'+str(ii)]=data['M059K_Rep'+str(ii)+'_T']
all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

for name in ['A172','GBM63']:
    data = pd.read_table('cell_lines/raw_data/'+name+'_U_CONvGAB_fpkm.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_CONGAB_'+str(ii)]=data['U_C_'+str(ii)]
        recurrent[name+'_CONGAB_'+str(ii)]=data['U_I_'+str(ii)]
    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

    data = pd.read_table('cell_lines/raw_data/'+name+'_U_CONvPTZ_fpkm.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_CONPTZ_'+str(ii)]=data['U_C_'+str(ii)]
        recurrent[name+'_CONPTZ_'+str(ii)]=data['U_I_'+str(ii)]
    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

    data = pd.read_table('cell_lines/raw_data/'+name+'_GAB_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_GAB_'+str(ii)]=data['U_'+str(ii)]
        recurrent[name+'_GAB_'+str(ii)]=data['T_'+str(ii)]
    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

    data = pd.read_table('cell_lines/raw_data/'+name+'_PTZ_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_PTZ_'+str(ii)]=data['U_'+str(ii)]
        recurrent[name+'_PTZ_'+str(ii)]=data['T_'+str(ii)]
    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

    data = pd.read_table('cell_lines/raw_data/'+name+'_CON_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_'+str(ii)]=data[name+'_U_'+str(ii)]
        recurrent[name+'_'+str(ii)]=data[name+'_T_'+str(ii)]
    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)
    
data = pd.read_table('cell_lines/raw_data/GBM63_CUTRUN2021_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
data=data[data.index.str.contains('_PAR_Y')==False]
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for ii in [0,1]:
    primary['GBM63_CUTRUN'+'_'+str(ii)]=data['U_'+str(ii)]
    recurrent['GBM63_CUTRUN'+'_'+str(ii)]=data['T_'+str(ii)]
all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week1_FPKM.txt',sep='\t',header=0,index_col=0)
data=data[data.index.str.contains('_PAR_Y')==False]
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for i in ['1','3']:
    data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week'+i+'_FPKM.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary['PDspheroids_'+i]=data.loc[:, data.columns.str.startswith('U')].mean(axis=1)
    recurrent['PDspheroids_'+i]=data.loc[:, data.columns.str.startswith('T')].mean(axis=1)
all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
data=data[data.index.str.contains('_PAR_Y')==False]
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for i in ['48','72']:
    data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
    data=data[data.index.str.contains('_PAR_Y')==False]
    data.index=[i[:15] for i in data.index]
    primary['ME'+'_'+i]=data.loc[:, data.columns.str.startswith('ME'+i+'hrU')].mean(axis=1)
    recurrent['ME'+'_'+i]=data.loc[:, data.columns.str.startswith('ME'+i+'hrT')].mean(axis=1)
print(primary.head())
print(recurrent.head())

all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

filt = pd.read_table('cell_lines/filtered_genes.txt',sep='\t',header=None,index_col=0)
all_prim=pd.merge(all_prim,filt, how = "inner",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,filt, how = "inner",left_index=True, right_index=True)

log2fc=pd.DataFrame(index=all_prim.index)

for c in all_prim.columns:
    prim=all_prim[c]
    recu=all_recu[c]
    log2fc[c]=list(np.log2((recu+0.01)/(prim+0.01)))
    log2fc[c].to_csv('cell_lines/ranks/actual_log2fc/'+c+'.rnk',sep='\t',header=False,index=True)
    abs(log2fc[c]).to_csv('cell_lines/ranks/absolute_log2fc/'+c+'.rnk',sep='\t',header=False,index=True)

log2fc.to_csv('cell_lines/tables/log2fc_cell_lines.txt',sep='\t',header=True,index=True)
all_prim.to_csv('cell_lines/tables/primary_cell_lines.txt',sep='\t',header=True,index=True)
all_recu.to_csv('cell_lines/tables/recurrent_cell_lines.txt',sep='\t',header=True,index=True)


