# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:14:44 2020

@author: medgnt
"""

import argparse
from os import listdir
from os import path

parser = argparse.ArgumentParser(description="Combine GSEA patient results.")
parser.add_argument('-g', '--genes', dest='genes', help='File with a list of genes, one per line.')
parser.add_argument('-d', '--directory', dest='directory', help='Directory containing results tables for each patient.')
parser.add_argument('-s', '--set', dest='set', help='Set.')

args = parser.parse_args()

nes={}
pval={}
fdr={}
patients=[]

with open(args.genes,'r') as file:
    for line in file:
        gene=line.strip()
        nes[gene]=[]
        pval[gene]=[]
        fdr[gene]=[]

for f in listdir(args.directory+'/'+args.set):
    with open(args.directory+'/'+args.set+'/'+f,'r') as file:
        p_nes={}
        p_pval={}
        p_fdr={}
        patient=path.basename(f)[:f.find('_')]
        patients.append(patient)
        for line in file:
            l=line.split('\t')
            if len(l)>1:
                p_nes[l[0]]=l[5]
                p_pval[l[0]]=l[6]
                p_fdr[l[0]]=l[7]
    for g in nes.keys():
        if g in p_nes:
            nes[g].append(p_nes[g])
            pval[g].append(p_pval[g])
            fdr[g].append(p_fdr[g])
        else:
            nes[g].append('NA')
            pval[g].append('NA')
            fdr[g].append('NA')

with open(args.directory+'/'+args.set+'_nes.txt','w+') as file:
    file.write('Gene\t'+'\t'.join(patients))
    for g in nes:
        file.write('\n'+g+'\t'+'\t'.join(nes[g]))
with open(args.directory+'/'+args.set+'_pval.txt','w+') as file:
    file.write('Gene\t'+'\t'.join(patients))
    for g in pval:
        file.write('\n'+g+'\t'+'\t'.join(pval[g]))
with open(args.directory+'/'+args.set+'_fdr.txt','w+') as file:
    file.write('Gene\t'+'\t'.join(patients))
    for g in fdr:
        file.write('\n'+g+'\t'+'\t'.join(fdr[g]))

