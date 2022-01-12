import pandas as pd
import math
import matplotlib.pyplot as plt

up={}
down={}
with open('deseq2_uvd/results_up.txt','r') as file1:
    file1.readline()
    for line in file1:
        l=line.strip().split(' ')
        if l[2]=='NA':
            sign=1
        elif float(l[2])<0:
            sign=-1
        else:
            sign=1
        if l[5]=='NA':
            val=1
        else:
            val=float(l[5])
        up[l[0].strip('"')[:15]]=sign*-math.log10(val)

with open('deseq2_uvd/results_down.txt','r') as file2:
    file2.readline()
    for line in file2:
        l=line.strip().split(' ')
        if l[2]=='NA':
            sign=1
        elif float(l[2])<0:
            sign=-1
        else:
            sign=1
        if l[5]=='NA':
            val=1
        else:
            val=float(l[5])
        down[l[0].strip('"')[:15]]=sign*-math.log10(val)


jbs=pd.read_table('gene_lists/JARID2_bound_genes.txt',sep='\t', header=None)
jbs=list(jbs[0].str.strip())
le50=pd.read_table('analysis/leading_edge/le50_actual_1000_gbm_idhwt_rt_tmz_local.txt',sep='\t', header=None)
le50=list(le50[0])
le70=pd.read_table('analysis/leading_edge/le70_actual_1000_gbm_idhwt_rt_tmz_local.txt',sep='\t', header=None)
le70=list(le70[0])

ks=[i for i in up.keys() if i not in jbs]
y=[up[i] for i in ks]
x=[down[i] for i in ks]
y1=[up[i[:15]] for i in jbs if i not in le50]
x1=[down[i[:15]] for i in jbs if i not in le50]
y2=[up[i] for i in le50 if i not in le70]
x2=[down[i] for i in le50 if i not in le70]
y3=[up[i] for i in le70]
x3=[down[i] for i in le70]

plt.figure()
plt.axhline(0,color='black',linestyle='dashed', linewidth=1,zorder=0)
plt.axvline(0,color='black',linestyle='dashed', linewidth=1,zorder=0)
plt.scatter(x,y, c='grey',edgecolors='none',s=20,alpha=0.6, label='All genes',zorder=10)
plt.scatter(x1,y1, c='red',edgecolors='none',s=20,alpha=0.6, label='JBS genes',zorder=10)
plt.scatter(x2,y2, c='blue',edgecolors='none',s=20,alpha=0.6, label='LE50 genes',zorder=10)
plt.scatter(x3,y3, c='limegreen',edgecolors='none',s=20,alpha=0.6, label='LE70 genes',zorder=10)
plt.legend()
plt.xlabel('Down DEA pvalues')
plt.ylabel('Up DEA pvalues')
plt.show()

plt.savefig('deseq2_uvd/pval_plot_LE.pdf')

pc=pd.read_table('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_head1000.txt',sep='\t', header=None)
pc=list(pc[0])
pc2=pd.read_table('analysis/pca/pca_gbm_idhwt_rt_tmz_local_log2fc_all_FALSE_PC1_loadings_head100.txt',sep='\t', header=None)
pc2=list(pc2[0])

ks=[i for i in up.keys() if i not in pc]
y=[up[i] for i in ks]
x=[down[i] for i in ks]
y1=[up[i[:15]] for i in pc if i not in pc2]
x1=[down[i[:15]] for i in pc if i not in pc2]
y2=[up[i[:15]] for i in pc2]
x2=[down[i[:15]] for i in pc2]

plt.figure()
plt.axhline(0,color='black',linestyle='dashed', linewidth=1,zorder=0)
plt.axvline(0,color='black',linestyle='dashed', linewidth=1,zorder=0)
plt.scatter(x,y, c='grey',edgecolors='none',s=20,alpha=0.6, label='All genes',zorder=10)
plt.scatter(x1,y1, c='red',edgecolors='none',s=20,alpha=0.6, label='PC1 top 1000 genes',zorder=10)
plt.scatter(x2,y2, c='blue',edgecolors='none',s=20,alpha=0.6, label='PC1 top 100 genes',zorder=10)
plt.legend()
plt.xlabel('Down DEA pvalues')
plt.ylabel('Up DEA pvalues')


plt.savefig('deseq2_uvd/pval_plot_PC1.pdf')

