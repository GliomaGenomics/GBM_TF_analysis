import pandas as pd
import numpy as np

all_prim=pd.DataFrame()
all_recu=pd.DataFrame()

for name in ['A172','GBM63']:
    data = pd.read_table('cell_lines/raw_data/'+name+'_CON_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    primary=pd.DataFrame(index=data.index)
    recurrent=pd.DataFrame(index=data.index)
    for ii in [0,1,2]:
        primary[name+'_'+str(ii)]=data[name+'_U_'+str(ii)]
        recurrent[name+'_'+str(ii)]=data[name+'_T_'+str(ii)]

    all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
    all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)
    
data = pd.read_table('cell_lines/raw_data/GBM63_CUTRUN2021_UvT_fpkm.txt',sep='\t',header=0,index_col=0)
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for ii in [0,1]:
    primary['GBM63_CUTRUN'+'_'+str(ii)]=data['U_'+str(ii)]
    recurrent['GBM63_CUTRUN'+'_'+str(ii)]=data['T_'+str(ii)]

all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week1_FPKM.txt',sep='\t',header=0,index_col=0)
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for i in ['1','3']:
    data = pd.read_table('cell_lines/raw_data/PDspheroids_singlecells_week'+i+'_FPKM.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    primary['PDspheroids_'+i]=data.loc[:, data.columns.str.startswith('U')].mean(axis=1)
    recurrent['PDspheroids_'+i]=data.loc[:, data.columns.str.startswith('T')].mean(axis=1)

all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
data.index=[i[:15] for i in data.index]
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)
for i in ['48','72']:
    data = pd.read_table('cell_lines/raw_data/ME_Oct2019_FPKM_all.txt',sep='\t',header=0,index_col=0)
    data.index=[i[:15] for i in data.index]
    primary['ME'+'_'+i]=data.loc[:, data.columns.str.startswith('ME'+i+'hrU')].mean(axis=1)
    recurrent['ME'+'_'+i]=data.loc[:, data.columns.str.startswith('ME'+i+'hrT')].mean(axis=1)

all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)

all_prim=all_prim[~all_prim.index.duplicated(keep='first')]
all_recu=all_recu[~all_recu.index.duplicated(keep='first')]

all_prim_filt=pd.DataFrame(index=all_prim.index)
all_recu_filt=pd.DataFrame(index=all_prim.index)
all_log2fc_filt=pd.DataFrame(index=all_prim.index)

for c in all_prim.columns:
    prim=pd.DataFrame(all_prim[c],index=all_prim.index)
    recu=pd.DataFrame(all_recu[c],index=all_recu.index)
    lqp=prim[c][prim[c]>0].quantile(0.25)
    lqr=recu[c][recu[c]>0].quantile(0.25)
    t=max([lqp,lqr])
    filt=((prim[c]>t)|(recu[c]>t)) 
#& (prim[c]>0) & (recu[c]>0)
    prim=prim[filt] 
    recu=recu[filt]
    vals=list(np.log2((recu[c]+0.01)/(prim[c]+0.01)))
    log2fc=pd.DataFrame(vals,columns=[c],index=prim.index)
    log2fc.to_csv('cell_lines/ranks/actual_log2fc/'+c+'.rnk',sep='\t',header=False,index=True)
    abs(log2fc).to_csv('cell_lines/ranks/absolute_log2fc/'+c+'.rnk',sep='\t',header=False,index=True)
    all_prim_filt=pd.merge(all_prim_filt,prim,how="outer",left_index=True, right_index=True).fillna('NA')
    all_recu_filt=pd.merge(all_recu_filt,recu,how="outer",left_index=True, right_index=True).fillna('NA')
    all_log2fc_filt=pd.merge(all_log2fc_filt,log2fc,how="outer",left_index=True, right_index=True).fillna('NA')

all_log2fc_filt.to_csv('cell_lines/tables/log2fc_cell_lines.txt',sep='\t',header=True,index=True)
all_prim_filt.to_csv('cell_lines/tables/primary_cell_lines.txt',sep='\t',header=True,index=True)
all_recu_filt.to_csv('cell_lines/tables/recurrent_cell_lines.txt',sep='\t',header=True,index=True)



