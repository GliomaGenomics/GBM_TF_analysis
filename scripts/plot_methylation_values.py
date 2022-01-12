import json
import seaborn as sns
import matplotlib.pyplot as plt

with open("methylation/methylation_values.json","r") as infile:
    vals=json.load(infile)


name={}
name['all_prim']='All genes in P'
name['all_recu']='All genes in R'
name['jar_prim']='JBS genes in P'
name['jar_recu']='JBS genes in R'
name['le50_prim']='LE50 genes in P'
name['le50_recu']='LE50 genes in R'
name['le70_prim']='LE70 genes in P'
name['le70_recu']='LE70 genes in R'

for i in vals:
    print([i,len(vals[i])])

line=["-","--","-","--","-","--","-","--"]
col=['black','black','darkgrey','darkgrey','#81bbff','#81bbff','#2943ff','#2943ff']
plt.figure()
i=0
for g in ['all_prim','all_recu','jar_prim','jar_recu','le50_prim','le50_recu','le70_prim','le70_recu']:
    sns.distplot(vals[g], hist = False, kde = True,kde_kws = {'linewidth': 2,"linestyle":line[i]},label = name[g], color=col[i])
    i+=1

plt.legend()
plt.xlabel('Promotor DNA methylation')
plt.ylabel('Frequency')
plt.savefig('methylation/methylation_plot.pdf')
