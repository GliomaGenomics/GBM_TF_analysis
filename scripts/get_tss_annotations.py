import pandas as pd

tss_transcript={}
tss_gene={}

with open('intermediate_files/transcript_to_tss_position.txt','r') as file:
    for line in file:
        l=line.strip().split('\t')
        if l[1] in tss_transcript:
            tss_transcript[l[1]].append(l[0])
        else:
            tss_transcript[l[1]]=[l[0]]
with open('intermediate_files/gene_to_tss_position.txt','r') as file:
    for line in file:
        l=line.strip().split('\t')
        if l[1] in tss_gene:
            tss_gene[l[1]].append(l[0].split('.')[0])
            print(l)
        else:
            tss_gene[l[1]]=[l[0].split('.')[0]]

print(stop)

mrna=list(pd.read_table('ranks/filtered_genelist_mrna.txt',sep='\t', header=None,index_col=None)[0])
total=list(pd.read_table('ranks/filtered_genelist_total.txt',sep='\t', header=None,index_col=None)[0])
jbs=list(pd.read_table('gene_lists/JARID2_bound_genes.txt',sep='\t', header=None,index_col=None)[0])
le50=list(pd.read_table('analysis/leading_edge/le50_actual_1000_gbm_idhwt_rt_tmz_local.txt',sep='\t', header=None,index_col=None)[0])
le70=list(pd.read_table('analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt',sep='\t', header=None,index_col=None)[0])
with open('gsea_files/TSS_TFs_ENS_1000_GTRDv19_10_gencodev27_JARID2_only.gmt','r') as file:
    jbstss=file.readline().strip().split('\t')[2:]  


with open('analysis/leading_edge/tss_annotations.txt','w+') as file:
    file.write('\t'.join(['tss','gene','transcripts','gene_filter','gene_status','jbs_tss'])+'\n')
    for t in tss_transcript.keys():
        for g in tss_gene[t]:
            out=[t,g,str(tss_transcript[t])]
            if g in mrna and g in total:
                out.append('all')
            elif g in mrna:
                out.append('mrna')
            elif g in total:
                out.append('total')
            else:
                out.append('none')
            if g in le70:
                out.append('le70')
            elif g in le50:
                out.append('le50')
            elif g in jbs:
                out.append('jbs')
            else:
                out.append('none')
            if t in jbstss:
                out.append('1')
            else:
                out.append('0')
            file.write('\t'.join(out)+'\n')
