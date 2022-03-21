import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Merge files on first column.")
parser.add_argument('-i', '--input', dest='input', help='Input.')
parser.add_argument('-o', '--output', dest='output', help='Output.')
args = parser.parse_args()

data=pd.read_table(args.input,sep='\t',header=0,index_col=0)

ens=pd.read_table('downloaded_data/gencode.v27_geneidtoname.txt',sep='\t', header=None)

#read in dictionary of gene to ensids
ensd={}
for r, en in ens.iterrows():
    ensd[en[0][0:15]]=en[1]

#convert data indexes from genenames to ensids
data=data.rename(ensd,axis='rows')

data.to_csv(args.output,sep='\t',header=True,index=True)
