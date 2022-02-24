import pandas as pd
import numpy as np

data = pd.read_table('tissue_slices/downloaded_data/GSE179649_Neuro_Organo_kallisto_gene_tpm.csv',sep=',',header=0,index_col=0)
data=data[~data.index.str.contains("_PAR_")]
data.index=[i[:15] for i in data.index]
patients=set([i[:7] for i in data.columns])

all_values=[]
for col in data:
    all_values.extend(list(data[col]))
all_values=pd.DataFrame(all_values)
lq=all_values[all_values >0].quantile(0.25)
print("lq="+str(lq))
ps=[]
for i in range(0,len(data)):
    abovep=data.iloc[i][data.iloc[i] >= lq[0]][~data.iloc[i][data.iloc[i] >= lq[0]].index.str.contains('TMZ')].count()
    totalp=data.iloc[i][~data.iloc[i].index.str.contains('TMZ')].count()
    proportionp=abovep/totalp
    abover=data.iloc[i][data.iloc[i] >= lq[0]][data.iloc[i][data.iloc[i] >= lq[0]].index.str.contains('TMZ')].count()
    totalr=data.iloc[i][data.iloc[i].index.str.contains('TMZ')].count()
    proportionr=abover/totalr
    if proportionp >=0.2 or proportionr >=0.2:
        ps.append('True_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))
    else:
        ps.append('False_'+str(proportionp)+'_'+str(abovep)+'_'+str(totalp)+'_'+str(proportionr)+'_'+str(abover)+'_'+str(totalr))

psd=pd.DataFrame(ps,columns=['ps'])
data=data[list(psd['ps'].str.contains('True'))]
log2fc=pd.DataFrame(index=data.index)
primary=pd.DataFrame(index=data.index)
recurrent=pd.DataFrame(index=data.index)

for p in patients:
    treated=[]
    untreated=[]
    for c in data.columns:
        if p in c:
            if 'TMZ_4Gy' in c:
                treated.append(c)
            else:
                untreated.append(c)
    for u in untreated:
        for t in treated:
            prim=data[u]
            recu=data[t]
            log2fc[str(u)+'_'+str(t)]=list(np.log2((recu+0.01)/(prim+0.01)))
            primary[str(u)+'_'+str(t)]=list(prim)
            recurrent[str(u)+'_'+str(t)]=list(recu)
            log2fc[str(u)+'_'+str(t)].to_csv('tissue_slices/ranks/actual_log2fc/'+str(u)+'_'+str(t)+'.rnk',sep='\t',header=False,index=True)
            abs(log2fc[str(u)+'_'+str(t)]).to_csv('tissue_slices/ranks/absolute_log2fc/'+str(u)+'_'+str(t)+'.rnk',sep='\t',header=False,index=True)
#   get summed samples across replicates if more than 1 u or t
    if len(untreated)+len(treated)>2:
        prim=data[untreated].sum(axis=1)
        recu=data[treated].sum(axis=1)
        log2fc[p+'_sum']=list(np.log2((recu+0.01)/(prim+0.01)))
        primary[p+'_sum']=list(prim)
        recurrent[p+'_sum']=list(recu)
        log2fc[p+'_sum'].to_csv('tissue_slices/ranks/actual_log2fc/'+p+'_sum.rnk',sep='\t',header=False,index=True)
        abs(log2fc[p+'_sum']).to_csv('tissue_slices/ranks/absolute_log2fc/'+p+'_sum.rnk',sep='\t',header=False,index=True)
log2fc.to_csv('tissue_slices/tables/log2fc_tissue_slices.txt',sep='\t',header=True,index=True)
primary.to_csv('tissue_slices/tables/primary_tissue_slices.txt',sep='\t',header=True,index=True)
recurrent.to_csv('tissue_slices/tables/recurrent_tissue_slices.txt',sep='\t',header=True,index=True)
