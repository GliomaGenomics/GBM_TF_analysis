import pandas
import seaborn as sns
import matplotlib.pyplot as plt

patients=list(pandas.read_table('patient_lists/glass_gbm_idhwt_rt_tmz_local_methylation.txt',header=None)[0])
beta=pandas.read_table('original_data/beta.merged.tsv', header=0)

all_pro_probes=[]
all_dict={}
with open('methylation/all_genes_probes.txt', 'r') as all_genes:
    for line in all_genes:
        l=line.strip().split()
        if len(l)>1:
            all_dict[l[0]]=l[1:]
            all_pro_probes.extend(l[1:])

jar_dict={}
with open('methylation/all_JARID2_bound_genes_probes.txt', 'r') as jar_genes:
    for	line in	jar_genes:
       	l=line.strip().split()
       	if len(l)>1:
       	    jar_dict[l[0]]=l[1:]

le50_dict={}
with open('methylation/le50_JARID2_bound_genes_probes.txt', 'r') as le50_genes:
    for line in le50_genes:
        l=line.strip().split()
        if len(l)>1:
            le50_dict[l[0]]=l[1:]

filt_beta=pandas.DataFrame()
filt_beta['probeID']=beta['probeID']
for p in patients:
    filt_beta[p+'-TP']=beta[beta.columns[beta.columns.str.startswith(p+'-TP')][0]]
    filt_beta[p+'-R1']=beta[beta.columns[beta.columns.str.startswith(p+'-R1')][0]]

print(beta.shape)
print(filt_beta.shape)

probes=pandas.DataFrame({'probeID':all_pro_probes})
filt_beta=filt_beta.merge(probes, how='inner', on='probeID')

print(filt_beta.shape)

del beta

vals={}
vals['all_prim']=[]
vals['all_recu']=[]
vals['jar_prim']=[]
vals['jar_recu']=[]
vals['le50_prim']=[]
vals['le50_recu']=[]
name={}
name['all_prim']='All genes in P'
name['all_recu']='All genes in R'
name['jar_prim']='JBS genes in P'
name['jar_recu']='JBS genes in R'
name['le50_prim']='LE50 genes in P'
name['le50_recu']='LE50 genes in R'

for gene in all_dict:
    probes=pandas.DataFrame({'probeID':all_dict[gene]})    
    fp=filt_beta.merge(probes, how='inner', on='probeID')
    if fp.shape[0]>0:
        for p in patients:
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-TP')][0]])
            vals['all_prim'].extend([sum(fpvals)/len(fpvals)])
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-R1')][0]])
            vals['all_recu'].extend([sum(fpvals)/len(fpvals)])

for gene in jar_dict:
    probes=pandas.DataFrame({'probeID':jar_dict[gene]})
    fp=filt_beta.merge(probes, how='inner', on='probeID')
    if fp.shape[0]>0:
        for p in patients:
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-TP')][0]])
            vals['jar_prim'].extend([sum(fpvals)/len(fpvals)])
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-R1')][0]])
            vals['jar_recu'].extend([sum(fpvals)/len(fpvals)])

for gene in le50_dict:
    probes=pandas.DataFrame({'probeID':le50_dict[gene]})
    fp=filt_beta.merge(probes, how='inner', on='probeID')
    if fp.shape[0]>0:
        for p in patients:
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-TP')][0]])
            vals['le50_prim'].extend([sum(fpvals)/len(fpvals)])
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-R1')][0]])
            vals['le50_recu'].extend([sum(fpvals)/len(fpvals)])

import json
json = json.dumps(vals)
f = open("methylation/methylation_values.json","w")
f.write(json)
f.close()

for i in vals:
    print([i,len(vals[i])])

line=["-","--","-","--","-","--"]
col=['black','black','#2943ff','#2943ff','#81bbff','#81bbff']
plt.figure()
i=0
for g in ['all_prim','all_recu','jar_prim','jar_recu','le50_prim','le50_recu']:
    sns.distplot(vals[g], hist = False, kde = True,kde_kws = {'linewidth': 2,"linestyle":line[i]},label = name[g], color=col[i])
    i+=1

plt.legend()
plt.xlabel('Promotor DNA methylation')
plt.ylabel('Frequency')
plt.savefig('methylation/methylation_plot.pdf')
