import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)

adata = ad.io.read_h5ad("/w5home/bmoore/Pierre_sc_zebrafish/output_seurat2ann_20260224_121120/seuratobj_all_reclustered_075.h5ad")
print(set(adata.obs['Type']))
target_types = ['1Wk', '2Wk-R', '4Wk-R', '6Wk-R']
rdata = adata[adata.obs['Type'].isin(target_types)].copy()
target_types2 = ['1Wk', '2Wk-NR', '4Wk-NR', '6Wk-NR']
nrdata = adata[adata.obs['Type'].isin(target_types2)].copy()

rdata.X = csr_matrix(rdata.X)
nrdata.X = csr_matrix(nrdata.X)

rdata.write_h5ad("rdata.h5ad")
nrdata.write_h5ad("nrdata.h5ad")