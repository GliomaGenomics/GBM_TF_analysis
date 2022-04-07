import pandas as pd

vars = pd.read_table('glass_data/variants_passgeno_20201109_filtered.csv',sep=',',header=None,index_col=None)
anos=pd.read_table('glass_data/variants_anno_20201109_silent.csv',sep=',',header=None,index_col=None)
anof=pd.read_table('glass_data/variants_anno_20201109_deleterious.csv',sep=',',header=None,index_col=None)
vars['id']=vars[1].astype(str)+':'+vars[2].astype(str)+'-'+vars[3].astype(str)+':'+vars[4].astype(str)
anos['id']=anos[0].astype(str)+':'+anos[1].astype(str)+'-'+anos[2].astype(str)+':'+anos[4].astype(str)
anof['id']=anof[0].astype(str)+':'+anof[1].astype(str)+'-'+anof[2].astype(str)+':'+anof[4].astype(str)
datas=vars.merge(anos, how='left', on='id')
dataf=vars.merge(anof, how='left', on='id')
datas=datas[datas['5_y'].notnull()]
dataf=dataf[dataf['5_y'].notnull()]
datas['sample']=datas['0_x'].str[:15]
dataf['sample']=dataf['0_x'].str[:15]
datas=datas[['sample','5_y']]
dataf=dataf[['sample','5_y']]
datas=datas.drop_duplicates()
dataf=dataf.drop_duplicates()

genes=pd.concat([datas['5_y'],dataf['5_y']]).unique()
 

meta=pd.read_table('reports/jarid2_results/outputs_actual_glass_1000_JARID2_results+dis.tsv', header=None, names = ['patient_id','value'],index_col=0)
patients=pd.read_table('patient_lists/glass_gbm_idhwt_rt_tmz_local+dis_dna.txt', header=None, index_col=None)

prim_up=[]
recu_up=[]
prim_down=[]
recu_down=[]

for p in patients[0]:
    if meta.at[p.replace('-','.'),'value']>0:
        prim_up.append(p+'-TP')
        recu_up.append(p+'-R1')
    elif meta.at[p.replace('-','.'),'value']<0:
        prim_down.append(p+'-TP')
        recu_down.append(p+'-R1')

groups=[prim_up,recu_up,prim_down,recu_down]
genegroups=[]
#append gene name for every variant in order: 
#prim_up_silent,prim_up_deleterious,recu_up_silent,recu_up_deleterious,prim_down_silent,prim_down_deleterious,recu_down_silent,recu_down_deleterious
for i in [0,1,2,3]:
    genegroups.append(list(datas[datas['sample'].isin(groups[i])]['5_y']))
    genegroups.append(list(dataf[dataf['sample'].isin(groups[i])]['5_y']))


with open('variants/variant_counts.txt','w+') as file:
    file.write('\t'.join(['gene','prim_up_silent','prim_up_deleterious','recu_up_silent','recu_up_deleterious','prim_down_silent','prim_down_deleterious','recu_down_silent','recu_down_deleterious'])+'\n')
    for g in genes:
        if g!='nan':
            l=[g]
            for c in genegroups:
                l.append(str(c.count(g)))
            file.write('\t'.join(l)+'\n')
