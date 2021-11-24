
import pandas as pd
import math
import matplotlib.pyplot as plt

up={}
down={}
upsign={}
downsign={}
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
        up[l[0].strip('"')[:15]]=-math.log10(val)
        upsign[l[0].strip('"')[:15]]=sign*-math.log10(val)

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
        down[l[0].strip('"')[:15]]=-math.log10(val)
        downsign[l[0].strip('"')[:15]]=sign*-math.log10(val)


#read in dictionary of gene ids to hgnc
names=pd.read_table('downloaded_data/gencode.v27_geneidtohgnc.txt',sep='\t', header=None)
ensd={}
for r, en in names.iterrows():
    ensd[en[0][0:15]]=str(en[1])

nameupens={}
namedownens={}
nameupenssign={}
namedownenssign={}
nameup={}
namedown={}
nameupsign={}
namedownsign={}
for i in up.keys():
    if i not in nameupens:
        if i in ensd:
            nameup[ensd[i]]=[]
            namedown[ensd[i]]=[]
            nameupsign[ensd[i]]=[]
            namedownsign[ensd[i]]=[]
        nameupens[i]=[]
        namedownens[i]=[]
        nameupenssign[i]=[]
        namedownenssign[i]=[]
    if i in ensd:
        nameup[ensd[i]].append(up[i])
        namedown[ensd[i]].append(down[i])
        nameupsign[ensd[i]].append(upsign[i])
        namedownsign[ensd[i]].append(downsign[i])
    nameupens[i].append(up[i])
    namedownens[i].append(down[i])
    nameupenssign[i].append(upsign[i])
    namedownenssign[i].append(downsign[i])


for i in nameup:
    nameup[i]=sum(nameup[i])/len(nameup[i])
    namedown[i]=sum(namedown[i])/len(namedown[i])
    nameupsign[i]=sum(nameupsign[i])/len(nameupsign[i])
    namedownsign[i]=sum(namedownsign[i])/len(namedownsign[i])
for i in nameupens:
    nameupens[i]=sum(nameupens[i])/len(nameupens[i])
    namedownens[i]=sum(namedownens[i])/len(namedownens[i])
    nameupenssign[i]=sum(nameupenssign[i])/len(nameupenssign[i])
    namedownenssign[i]=sum(namedownenssign[i])/len(namedownenssign[i])

with open('deseq2_uvd/u_ranks_ens.rnk','w+') as file:
    for i in nameupens:
        file.write(i+'\t'+str(nameupens[i])+'\n')

with open('deseq2_uvd/d_ranks_ens.rnk','w+') as file:
    for i in nameupens:
        file.write(i+'\t'+str(namedownens[i])+'\n')

with open('deseq2_uvd/usign_ranks_ens.rnk','w+') as file:
    for i in nameupens:
        file.write(i+'\t'+str(nameupenssign[i])+'\n')

with open('deseq2_uvd/dsign_ranks_ens.rnk','w+') as file:
    for i in nameupens:
        file.write(i+'\t'+str(namedownenssign[i])+'\n')

with open('deseq2_uvd/u_ranks.rnk','w+') as file:
    for i in nameup:
        file.write(i+'\t'+str(nameup[i])+'\n')

with open('deseq2_uvd/d_ranks.rnk','w+') as file:
    for i in nameup:
        file.write(i+'\t'+str(namedown[i])+'\n')

with open('deseq2_uvd/usign_ranks.rnk','w+') as file:
    for i in nameup:
        file.write(i+'\t'+str(nameupsign[i])+'\n')

with open('deseq2_uvd/dsign_ranks.rnk','w+') as file:
    for i in nameup:
        file.write(i+'\t'+str(namedownsign[i])+'\n')


