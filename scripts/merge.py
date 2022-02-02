import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Merge files on first column.")
parser.add_argument('-i', '--input', dest='input', help='List of comma separated files to merge.')
parser.add_argument('-o', '--output', dest='output', help='Output.')
args = parser.parse_args()

data=pd.DataFrame()
for f in args.input.split(','):
    data2 = pd.read_table(f,sep='\t',header=0,index_col=0)
    data=data.merge(data2, how='outer', left_index=True, right_index=True)

data.index.name='Rows'
data.to_csv(args.output,sep='\t',header=True,index=True)
