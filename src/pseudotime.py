# Palantir
# source venv/bin/activate
# load packages
print("loading packages")

import palantir
import scanpy as sc
import pandas as pd
import os
import shutil
import numpy as np
import cellrank as cr
import scvelo as scv
import pandas as pd

# Plotting
import matplotlib
import matplotlib.pyplot as plt

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
#warnings.filterwarnings(
#    action="ignore", module="scanpy", message="No data for colormapping"
#)

# read config file
import json
GIT_DIR = os.getcwd()
with open(GIT_DIR+'/config.json') as f:
    config_dict = json.load(f)
    print("loaded config file: ", config_dict["pseudotime"])

# variables
DATA_DIR = config_dict["pseudotime"]["DATA_DIR"]
ADATA_FILE = config_dict["pseudotime"]["ADATA_FILE"] 
NC = config_dict["pseudotime"]["NC"] # number of components
METADATA = config_dict["pseudotime"]["METADATA"]
annot_label = config_dict["pseudotime"]["annot_label"]

# load data
## note- data was previously converted from seurat object to anndata object
data_dir = os.path.expanduser(DATA_DIR)
print("loading data")
adata = sc.read_h5ad(data_dir + ADATA_FILE)
print(adata)
from datetime import datetime
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")
out_dir = data_dir + "pseudotime_" + now +"/"
os.mkdir(out_dir)
# copy config file
shutil.copy(GIT_DIR+'/config.json', DATA_DIR) 

# add in metadata
if METADATA != "NA":
    metadata = pd.read_csv(DATA_DIR + METADATA, sep="\t", index_col=0)
    meta = metadata[annot_label]
    adata.obs[annot_label] = meta
else:
    print("no metadata loaded")
# check metadata
print("check metadata")
print(adata.obs.columns)
print(adata.obs[annot_label].unique())

# visualize (note we have already done nearest neighbor and umap)
with plt.rc_context():
    sc.pl.embedding(
        adata,
        basis="umap",
        frameon=False,
        color=annot_label,
        )
    plt.savefig(out_dir + "labeled_clusters.pdf")

## gene expression visualization on UMAP
print("visualize gene expression")
# get inputs from config
genes = config_dict["pseudotime"]["gene_list"]
magic_impute = config_dict["pseudotime"]["magic_impute"]

# MAGIC imputation
# Palantir can use MAGIC to impute the data for visualization and determining gene expression trends.
if magic_impute == "TRUE":
    print("Running MAGIC imputation")
    imputed_X = palantir.utils.run_magic_imputation(adata)
    # umap with magic
    with plt.rc_context():
        sc.pl.embedding(
        adata,
        basis="umap",
        layer="MAGIC_imputed_data",
        color=genes,
        frameon=False,
        show=False,
        )
    plt.savefig(out_dir + "gene_expression_magic.pdf")
else:
    print("no magic imputation, showing true gene expression only")

# umap with gene expression
with plt.rc_context():
    sc.pl.embedding(
    adata,
    basis="umap",
    use_raw=False,
    color=genes,
    frameon=False,
    show=False,
    )
    plt.savefig(out_dir + "gene_expression.pdf")

# Run diffusion maps
# do pca first just in case
sc.pp.pca(adata)
print("Running diffusion maps")
dm_res = palantir.utils.run_diffusion_maps(adata, n_components=NC)
ms_data = palantir.utils.determine_multiscale_space(adata)

# Diffusion maps visualization
palantir.plot.plot_diffusion_components(adata)
plt.savefig(out_dir + 'palantir_components_umap.pdf')

# find terminal states
print("Finding terminal cell types")
terminal_celltypes = config_dict["pseudotime"]["terminal_celltypes"]
terminal_states = palantir.utils.find_terminal_states(adata, celltypes=terminal_celltypes, celltype_column=annot_label)
print(terminal_states)

# find start cell
print("Finding start cell")
start_celltype = config_dict["pseudotime"]["start_celltype"]
start_cell = palantir.utils.early_cell(adata,celltype=start_celltype, celltype_column=annot_label, 
                                       fallback_seed=1234)

# use cells found for start cell and terminal states for palantir analysis

# Run Palantir
print("Running Palantir")
if len(terminal_states) == 0:
    lista = [start_cell]
    pr_res = palantir.core.run_palantir(adata, start_cell, num_waypoints=500, use_early_cell_as_start=True)
    # plot with start cell
    palantir.plot.highlight_cells_on_umap(adata, lista)
    plt.savefig(out_dir + 'palantir_start_cell.pdf')
else:
    pr_res = palantir.core.run_palantir(
    adata, start_cell, num_waypoints=500, terminal_states=terminal_states
    )
    # plot with terminal cells and start cell
    all_cells = terminal_states.copy()
    all_cells[start_cell] = "Start cell"
    palantir.plot.highlight_cells_on_umap(adata, all_cells)
    plt.savefig(out_dir + 'palantir_terminal_cells.pdf')
    print(all_cells)

# plot results
palantir.plot.plot_palantir_results(adata, s=3)
plt.savefig(out_dir + 'palantir_results.pdf')

# Terminal state probability distributions of individual cells can be visualized using the plot_terminal_state_probs function
print("Determining terminal state probabilities")
if len(terminal_states) == 0:
    lista = [start_cell]
    palantir.plot.plot_terminal_state_probs(adata, lista)
    plt.savefig(out_dir + 'terminal_state_probs.pdf')

else:
    cells=all_cells.index
    print(cells)
    palantir.plot.plot_terminal_state_probs(adata, cells)
    plt.tight_layout()
    plt.savefig(out_dir +'terminal_state_probs.pdf', bbox_inches='tight')

# Branch probabilities
print("Determining branch state probabilities")
# Before computing the gene expression trends, we first need to select cells associated with a specific branch of the pseudotime trajectory.
masks = palantir.presults.select_branch_cells(adata, eps=0)
# visualize the branch selection
palantir.plot.plot_branch_selection(adata)
plt.savefig(out_dir + 'branch_selection.pdf')

# Compute gene trends over pseudotime
print("Computing gene trends")
if magic_impute == "TRUE":
    gene_trends = palantir.presults.compute_gene_trends(
        adata, expression_key="MAGIC_imputed_data")#gene_trend_key
else:
    gene_trends = palantir.presults.compute_gene_trends(
        adata)
print(adata)

# plot gene trends
palantir.plot.plot_gene_trends(adata, genes)
plt.savefig(out_dir + 'gene_trends.pdf')

# Heatmap trend visualization
palantir.plot.plot_gene_trend_heatmaps(adata, genes)
plt.savefig(out_dir + 'gene_trends_heatmap.pdf')

# Clustering
print("Clustering gene trends for each terminal state")
# Gene expression trends can be clustered and visualized
# Taking variable genes
more_genes = adata.var_names
# cluster for each terminal cell type
for state in terminal_celltypes:
    try:
        communities = palantir.presults.cluster_gene_trends(adata, state, more_genes)
        palantir.plot.plot_gene_trend_clusters(adata, state)
        plt.savefig(out_dir + 'gene_trend_clusters_'+ state + '.pdf')
    except(KeyError):
        print(state, " not in terminal cell types")

# weird seurat to anndata thing- anndata doesn not like a column named _index
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

# compute diffusion pseudotime (DPT)
print("Computing nearest neighbors and diffusion pseudotime (DPT)")
# start by computing a diffusion map
# need to recompute nearest neighbors first
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
sc.tl.diffmap(adata)
# find root cell
root_ixs = adata.obsm['X_diffmap'][:, 3].argmax()
print(root_ixs)
# visualize and set root
scv.pl.scatter(
    adata,
    basis="diffmap",
    c=[annot_label, root_ixs],
    legend_loc="right",
    components=["2, 3"],
    save= out_dir +'diffusion_map.png'
    )
adata.uns["iroot"] = root_ixs
# Compute DPT and compare it with the precomputed Palantir pseudotime:
sc.tl.dpt(adata)
with plt.rc_context():
    sc.pl.embedding(
    adata,
    basis="umap",
    color=["dpt_pseudotime", "palantir_pseudotime"],
    color_map="gnuplot2",
    show=False
    )
    plt.savefig(out_dir +'dpt-palantir_pseudotime.pdf')
    
# plot trajectory for each terminal cell type
print("Plotting trajectories for each terminal cell type")
# get trajectory dict
trajectory_dict = config_dict["pseudotime"]["trajectory"]
for key in trajectory_dict:
    try:
        trajectory = trajectory_dict[key]
        mask = np.in1d(adata.obs[annot_label], trajectory)
        with plt.rc_context():
            sc.pl.violin(
            adata[mask],
            keys=["dpt_pseudotime", "palantir_pseudotime"],
            groupby=annot_label,
            rotation=-90,
            order=trajectory,
            show=False
            )
            plt.tight_layout()
            plt.savefig(out_dir + key + '_trajectory.pdf', bbox_inches='tight')
    except(KeyError):
        print(key, "not found, continuing")

# Compute a transition matrix based on Palantir pseudotime
print("Computing transition matrix based on Palantir pseudotime")
pk = cr.kernels.PseudotimeKernel(adata, time_key="palantir_pseudotime")
pk.compute_transition_matrix()
print(pk)
# visualize based on pseudotime and transition matrix
with plt.rc_context():
    pk.plot_projection(basis="umap", recompute=True, 
                       show=False)
    plt.savefig(out_dir +'pseudotime_transition.png')
    
# Save results
print("Saving results")
file_path = os.path.join(out_dir, "pseudotime_processed.h5ad")
adata.write(file_path)
# package versions
import pkg_resources
packagever = open(out_dir+"package_versions.txt", "w")
for package in pkg_resources.working_set:
    print(package.key, pkg_resources.get_distribution(package).version)
    packagever.write(package.key + " " + pkg_resources.get_distribution(package).version + "\n")
packagever.close()