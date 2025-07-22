#This script takes a list of genes and returns a list of TFs that are in the list
import sys, os
import pandas as pd

TF_list = sys.argv[1]
DF1 = sys.argv[2]
df1 = pd.read_csv(DF1, sep='\t')  # Do not set index_col, keep all columns
print(df1['human_symbol'].tolist())
TF_list = pd.read_csv(TF_list, sep='\t', index_col = 0)
print(list(TF_list.index))
TF_list = list(TF_list.index)
TF_list = [x.upper() for x in TF_list]

# Subset df1 where 'human_symbol' is in TF_list (case-insensitive)
df1['human_symbol_upper'] = df1['human_symbol'].str.upper()
df1_subset = df1[df1['human_symbol_upper'].isin(TF_list)]

print(df1_subset)
# Optionally, save to file:
df1_subset.to_csv('TFs_subset_zebrafish.tsv', sep='\t', index=False)




