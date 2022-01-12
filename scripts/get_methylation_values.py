import pandas

patients=list(pandas.read_table('patient_lists/glass_gbm_idhwt_rt_tmz_local_methylation+rna.txt',header=None)[0])
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
le70_dict={}
with open('methylation/le70_JARID2_bound_genes_probes.txt', 'r') as le70_genes:
    for line in le70_genes:
        l=line.strip().split()
        if len(l)>1:
            le70_dict[l[0]]=l[1:]

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
vals['le70_prim']=[]
vals['le70_recu']=[]

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
for gene in le70_dict:
    probes=pandas.DataFrame({'probeID':le70_dict[gene]})
    fp=filt_beta.merge(probes, how='inner', on='probeID')
    if fp.shape[0]>0:
        for p in patients:
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-TP')][0]])
            vals['le70_prim'].extend([sum(fpvals)/len(fpvals)])
            fpvals=list(fp[fp.columns[fp.columns.str.startswith(p+'-R1')][0]])
            vals['le70_recu'].extend([sum(fpvals)/len(fpvals)])

