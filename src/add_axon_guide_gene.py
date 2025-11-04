import sys, os
import pandas as pd
import numpy as np

gene_file = sys.argv[1]
DF1 = sys.argv[2]
df1 = pd.read_csv(DF1, sep='\t')  # Do not set index_col, keep all columns
print(df1.head())
gene_df = pd.read_csv(gene_file, sep='\t')
print(gene_df['Mouse gene name'].tolist())
gene_list = set(gene_df['Mouse gene name'].tolist())

#df1['axon_guide_gene'] = np.where(df1['ligand'].isin(gene_list) | df1['receptor'].isin(gene_list), True, False)
df1['secreted_protein'] = np.where(df1['ligand'].isin(gene_list) | df1['receptor'].isin(gene_list), "yes", "no")
df1.head()
DF1 =DF1.replace('.txt', '')
DF1 =DF1.replace('.tsv', '')
#df1.to_csv(DF1 + '_axonguides.tsv', sep='\t', index=False)
df1.to_csv(DF1 + '_secreted_prot.tsv', sep='\t', index=False)
