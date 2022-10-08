import pandas as pd
import numpy as np

all_prim=pd.DataFrame()
all_recu=pd.DataFrame()

filt = pd.read_table('ranks/actual_log2fc/AD20.rnk',sep='\t',header=None,index_col=0)
del filt[filt.columns[0]]


for name in ['A172','GBM63']:

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


all_prim=pd.merge(all_prim,primary, how = "outer",left_index=True, right_index=True)
all_recu=pd.merge(all_recu,recurrent, how = "outer",left_index=True, right_index=True)


all_values=[]
for col in all_prim:
    all_values.extend(list(all_prim[col]))
for col in all_recu:
    all_values.extend(list(all_recu[col]))
all_values=pd.DataFrame(all_values)
lq=all_values[all_values >0].quantile(0.25)
print("lq="+str(lq))
ps=[]
for i in range(0,len(all_prim)):
    abovep=all_prim.iloc[i][all_prim.iloc[i] >= lq[0]].count()
    totalp=all_prim.iloc[i].count()
    proportionp=abovep/totalp
    abover=all_recu.iloc[i][all_recu.iloc[i] >= lq[0]].count()
    totalr=all_recu.iloc[i].count()
    proportionr=abover/totalr
    if proportionp >=1 or proportionr >=1:
        ps.append('True_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
    else:
        ps.append('False_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))

psd=pd.DataFrame(ps,columns=['ps'])
all_prim=all_prim[list(psd['ps'].str.contains('True'))]
all_prim.to_csv('cell_lines/filtered_genes.txt', columns=[], header=False, sep='\t',index=True)



