# conda activate realtime

# import packages
print("loading packages")
import sys
#!pip install moscot
#import moscot
import matplotlib.pyplot as plt
# warnings
import warnings
import scanpy as sc, anndata as ad, numpy as np, pandas as pd
from scipy import sparse
from anndata import AnnData
from moscot.problems.time import TemporalProblem
import moscot.plotting as mtp
import cellrank as cr
from cellrank.kernels import RealTimeKernel
import os
import shutil
import moscot as mt
import anndata as ad
import h5py
import numpy as np
import subprocess
import sys

sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

# read config file
import json
#GIT_DIR = '/w5home/bmoore/scRNAseq_downstream/'
with open('config.json') as f:
    config_dict = json.load(f)
    print("loaded config file: ", config_dict["realtime"])
    
# variables
DATA_DIR = config_dict["realtime"]["DATA_DIR"]
ADATA_FILE = config_dict["realtime"]["ADATA_FILE"]
time_label = config_dict["realtime"]["time_label"]
annot_label = config_dict["realtime"]["annot_label"]
dim_red = config_dict["realtime"]["dim_red"]

# make out dir
data_dir = os.path.expanduser(DATA_DIR)
from datetime import datetime
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")
out_dir = data_dir + "realtime_" + now +"/"
os.mkdir(out_dir)
# copy config file
shutil.copy(GIT_DIR + '/config.json', out_dir) 

# load data
## note- data was previously converted from seurat object to anndata object
print("loading data")
if dim_red == "pca" or dim_red == "umap":
    adata = sc.read_h5ad(data_dir + ADATA_FILE)
# if adding a new dim_red component
elif dim_red != "NA" and dim_red != "":
    # load h5ad file
    adata = ad.read_h5ad(data_dir + ADATA_FILE)
    # replace .h5ad with .h5Seurat
    ADATA_FILE = ADATA_FILE.replace(".h5ad", ".h5Seurat")
    # Open the original h5Seurat file
    with h5py.File(data_dir + ADATA_FILE, "r") as f:
        # Extract the dim_red data and transpose it
        dim_red_component = f[dim_red+"1"][:].T
        # Add the integrated.cca component
    adata.obsm["X_new_dim_red"] = dim_red_component
    adata.write_h5ad(out_dir +"anndata_obj_with_"+ dim_red + ".h5ad")
else:
    adata = sc.read_h5ad(data_dir + ADATA_FILE)
print(adata)

adata.obs["time"] = adata.obs[time_label]
print(adata.obs["time"].unique())

# checking all elements to be numeric using isdigit()
res = all(ele.isdigit() for ele in adata.obs["time"].unique())

# replacing values to number if category
if res==True:
    pass
else:
    count=1
    for i in adata.obs["time"].unique():
        print(i)
        adata.obs["time"].replace([i],[count], inplace=True)
        count= count + 1
# convert to float
adata.obs["time"] = adata.obs["time"].astype(float).astype("category")
print("time points: ")
print(adata.obs["time"].unique())

# convert annotations to string
adata.obs[annot_label] = adata.obs[annot_label].astype("string").astype("category")
print("cell types or clusters:")
print(adata.obs[annot_label].unique())

# force directed layout 
print("Draw force-directed layout")
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30, random_state=0)
sc.tl.draw_graph(adata, layout='fa')
print(adata.obsm.keys())
try:
    with plt.rc_context():
        sc.pl.embedding(
            adata,
            basis="X_draw_graph_fa",
            color=["time", annot_label],
            color_map="gnuplot",
            show=False)
        plt.tight_layout()
        plt.savefig(out_dir + "force_directed_graph.pdf", bbox_inches='tight')
except KeyError:
    with plt.rc_context():
        sc.pl.embedding(
            adata,
            basis="X_draw_graph_fr",
            color=["time", annot_label],
            color_map="gnuplot",
            show=False)
        plt.tight_layout()
        plt.savefig(out_dir + "force_directed_graph.pdf", bbox_inches='tight')
# if data was merged, should redo variable genes and pca

# get proliferation and apopstosis genes
print("get proliferation and apopstosis genes")
from moscot.utils.data import proliferation_markers, apoptosis_markers
p_markers = proliferation_markers("human")
a_markers = apoptosis_markers("human")
print(p_markers,a_markers)
type(p_markers)

# get variable genes
print("get variable genes")
gene_pool = list(adata.var_names)
print(len(gene_pool))

# With moscot, we couple cells across time points using optimal transport (OT) by setting up the temporal problem.
print("set up optimal transport")
tp = TemporalProblem(adata)
# Next, we adjust the marginals for cellular growth- and death rates.
tp = tp.score_genes_for_marginals(
    gene_set_proliferation="human", gene_set_apoptosis="human",use_raw=False
)

# visualize the proliferation and apoptosis scores
with plt.rc_context():
    sc.pl.embedding(
        adata, basis="X_draw_graph_fa", color=["time", "proliferation", "apoptosis"],
        show=False)
    plt.tight_layout()
    plt.savefig(out_dir + "prolif-apop_graph.pdf", bbox_inches='tight')
    
# Following the original Waddington OT publication, 
# we use local PCAs, computed separately for each pair of time points, 
# to calculate distances among cells.
print("Calculate distances among cells for each time point pair")
if dim_red == "pca":
    tp = tp.prepare(time_key="time", joint_attr="X_pca")
elif dim_red == "umap":
    tp = tp.prepare(time_key="time", joint_attr="X_umap")
elif dim_red != "NA" and dim_red != "":
    tp = tp.prepare(time_key="time", joint_attr="X_new_dim_red")
else:
    tp = tp.prepare(time_key="time", joint_attr="X_pca")
# We solve one OT problem per time point pair, probabilistically matching 
# early to late cells.
print("Probalistically matching early to late cells")
tp = tp.solve(epsilon=1e-3, tau_a=0.95, tau_b=0.999, scale_cost="mean")
# Above, epsilon and tau_a control the amount of entropic regularization and unbalancedness on the source marginal, respectively. 
# Higher entropic regularization speeds up the optimization and improves statistical properties of the solution [Cuturi, 2013]; 
# unbalancedness makes the solution more robust with respect to uncertain cellular growth rates and biased cell sampling

# compare prior and posterior growth rates
print("compare prior and posterior growth rates")
adata.obs["prior_growth_rates"] = tp.prior_growth_rates
adata.obs["posterior_growth_rates"] = tp.posterior_growth_rates
with plt.rc_context():
    sc.pl.embedding(
        adata,
        basis="X_draw_graph_fa",
        color=["prior_growth_rates", "posterior_growth_rates"],
        vmax="p99",
    )
    plt.tight_layout()
    plt.savefig(out_dir + "prior-post_growthrates.pdf", bbox_inches='tight')
    
# add cell costs
print("add cell costs")
# High values indicate that a certain cell is unlikely to have a descendant or ancestor
adata.obs["cell_costs_source"] = tp.cell_costs_source
adata.obs["cell_costs_target"] = tp.cell_costs_target
# visulaize cell costs
with plt.rc_context():
    sc.pl.embedding(
        adata, basis="X_draw_graph_fa", 
        color=["cell_costs_source", "cell_costs_target"]
    )
    plt.tight_layout()
    plt.savefig(out_dir + "cell_costs.pdf", bbox_inches='tight')
    
# identify ancestry of cells

# investigate which ancestry population a certain cell type has. 
# We do this by aggregating the transport matrix by cell type, using cell_transition()
# forward=True means we are plotting descendents, 
# forward=False means we are plotting ancestors
# loop through timepoints
print("identify ancestry and descendants of cells")
timepoints = adata.obs["time"].unique()
print(timepoints)

for i in timepoints:
    for j in timepoints:
        if i < j:  # Only process if i is earlier than j
            print(i, j)
            key = f"transitions_{i}_{j}"
            tp.cell_transition(
                source=i, 
                target=j, 
                source_groups=annot_label, 
                target_groups=annot_label, 
                forward=True, 
                key_added=key
            )
            
            # Create plots
            mtp.cell_transition(adata, key=key, dpi=100, 
                fontsize=8, save=f"{out_dir}{key}.png")
            
            # Get all unique cell types across both timepoints
            all_celltypes = set(adata[adata.obs['time'].isin([i, j])].obs[annot_label].unique())
            
            for ct in all_celltypes:
                try:
                    # Ancestors
                    tp.pull(source=i, target=j, data=annot_label, subset=ct)
                    fig, axes = plt.subplots(ncols=2, figsize=(20, 6))
                    
                    # Check if the cell type exists in the target timepoint
                    target_cells = adata[(adata.obs['time'] == j) & (adata.obs[annot_label] == ct)]
                    if len(target_cells) > 0:
                        mtp.pull(
                            tp,
                            time_points=[j],
                            basis="X_draw_graph_fa",
                            ax=axes[0],
                            title=[f"{ct} at timepoint {j}"],
                        )
                    else:
                        axes[0].text(0.5, 0.5, f"No {ct} cells at timepoint {j}", 
                                     ha='center', va='center')
                        axes[0].set_title(f"{ct} at timepoint {j}")
                    
                    mtp.pull(
                        tp,
                        time_points=[i],
                        basis="X_draw_graph_fa",
                        ax=axes[1],
                        title=[f"{ct} ancestors"],
                    )
                    
                    fig.suptitle(f"Ancestors of {ct}")
                    plt.tight_layout()
                    plt.savefig(f"{out_dir}{ct}_ancestors_{j}-{i}.pdf", bbox_inches='tight')
                    plt.close(fig)
                    
                    # Descendants
                    tp.push(source=i, target=j, data=annot_label, subset=ct)
                    fig, axes = plt.subplots(ncols=2, figsize=(20, 6))
                    
                    # Check if the cell type exists in the source timepoint
                    source_cells = adata[(adata.obs['time'] == i) & (adata.obs[annot_label] == ct)]
                    if len(source_cells) > 0:
                        mtp.push(
                            tp,
                            time_points=[i],
                            basis="X_draw_graph_fa",
                            ax=axes[0],
                            title=[f"{ct} at time {i}"],
                        )
                    else:
                        axes[0].text(0.5, 0.5, f"No {ct} cells at timepoint {i}", 
                                     ha='center', va='center')
                        axes[0].set_title(f"{ct} at time {i}")
                    
                    mtp.push(
                        tp,
                        time_points=[j],
                        basis="X_draw_graph_fa",
                        ax=axes[1],
                        title=[f"{ct} descendants"],
                    )
                    
                    fig.suptitle(f"Descendants of {ct}")
                    plt.tight_layout()
                    plt.savefig(f"{out_dir}{ct}_descendants_{i}-{j}.pdf", bbox_inches='tight')
                    plt.close(fig)
                    
                except Exception as e:
                    print(f"Error processing {ct} between timepoints {i} and {j}: {str(e)}")
                    continue
            
            # Sankey diagram
            tp.sankey(
                source=i,
                target=j,
                source_groups=annot_label,
                target_groups=annot_label,
                threshold=0.05,
                forward=True,
            )
            mtp.sankey(tp, dpi=100, figsize=(8, 4), title=f"Cell type evolution {i} to {j}", 
                       save=f"{out_dir}cell_type_evolution{i}-{j}.pdf")

# get cell transition matrix
print(adata.uns['moscot_results']['cell_transition'])

# Set up the RealTimeKernel
print("compute transition matrix")
tmk = RealTimeKernel.from_moscot(tp)

# to get from OT transport maps to a markov chain:
# 1. we sparsify OT transport maps by removing entries below a certain threshold; entropic regularization yields dense matrices which would make CellRank analysis very slow.
# 2. we use OT transport maps and molecular similarity to model transitions across and within time points, respectively.
# 3. we row-normalize the resulting cell-cell transition matrix (including all time points) and construct the Markov chain.
tmk.compute_transition_matrix(self_transitions="all", conn_weight=0.2, threshold="auto")

# Visualize the recovered dynamics by sampling random walks.
# This method simulates random walks on the Markov chain defined though the corresponding transition matrix
# Random walks are simulated by iteratively choosing the next cell based on the current cell’s transition probabilities.
with plt.rc_context():
    tmk.plot_random_walks(
        max_iter=500,
        start_ixs={"time": 1},
        basis="X_draw_graph_fa",
        seed=0,
        dpi=150,
        size=30,
    )
    plt.tight_layout()
    plt.savefig(out_dir + "random_walks.pdf", bbox_inches='tight')
    
# Black and yellow dots denote random walks starting and finishing points, respectively.

# plot the probability mass flow in time
# Visualize outgoing flow from a cluster of cells
# first get cells only in first timepoint
firsttime_celltypes = set(adata[adata.obs['time'].isin([1.0])].obs[annot_label].unique())
# loop through to show prob mass flow
for ct in firsttime_celltypes:
    with plt.rc_context():
        plt.figure(figsize=(10,4))
        ax = tmk.plot_single_flow(
            cluster_key=annot_label,
            time_key="time",
            cluster=ct,
            min_flow=0.1,
            xticks_step_size=4,
            show=False,
            #clusters=celltype_list,
            )

    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.tight_layout()
    plt.savefig(out_dir + "prob_mass_flow_" + ct + ".pdf", bbox_inches='tight')
    
# save object
print("saving object and package versions")
import anndata as ad
import pandas as pd
import numpy as np

# convert annotation to numpy array in order to save
adata.obs[annot_label] = np.array(adata.obs[annot_label])
# Verify the written file
try:
    adata_read = ad.read_h5ad(out_dir + "realtime_object.h5ad")
    print("Successfully read the written AnnData object.")
except Exception as e:
    print(f"Error reading the written AnnData object: {str(e)}")
# now save
adata.write(out_dir + "realtime_object.h5ad")
# save package versions
# Check if we're in a conda environment
in_conda = os.environ.get('CONDA_DEFAULT_ENV') is not None
if in_conda:
    # Use conda list command
    result = subprocess.run(['conda', 'list', '--explicit'], capture_output=True, text=True)
    packages = result.stdout
    # Write the packages to a file
    with open(out_dir + 'conda-requirements.txt', 'w') as f:
        f.write(packages)
    print("Package versions have been written to conda-requirements.txt")
else:
    print("Not in a conda environment. Please activate your environment first.")