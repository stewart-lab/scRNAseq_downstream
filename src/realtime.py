# conda activate realtime

# import packages
print("loading packages")
import sys
#!pip install moscot
#import moscot
import matplotlib.pyplot as plt
# warnings
import warnings

from moscot.problems.time import TemporalProblem
import moscot.plotting as mtp
import cellrank as cr
import scanpy as sc
from cellrank.kernels import RealTimeKernel
import os
import shutil
import moscot as mt

sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2

# read config file
import json
GIT_DIR = '/w5home/bmoore/scRNAseq_downstream/'
with open(GIT_DIR+'config.json') as f:
    config_dict = json.load(f)
    print("loaded config file: ", config_dict["realtime"])
    
# variables
DATA_DIR = config_dict["realtime"]["DATA_DIR"]
ADATA_FILE = config_dict["realtime"]["ADATA_FILE"]
time_label = config_dict["realtime"]["time_label"]
annot_label = config_dict["realtime"]["annot_label"]

# load data
## note- data was previously converted from seurat object to anndata object
data_dir = os.path.expanduser(DATA_DIR)
print("loading data")
adata = sc.read_h5ad(data_dir + ADATA_FILE)
print(adata)

# make out dir
from datetime import datetime
now = datetime.now()
now = now.strftime("%Y%m%d_%H%M%S")
out_dir = data_dir + "realtime_" + now +"/"
os.mkdir(out_dir)
# copy config file
shutil.copy(GIT_DIR+'config.json', out_dir) 

adata.obs["time"] = adata.obs[time_label]#.astype(float).astype("category")
print(adata.obs["time"].unique())

# replacing values
count=1
for i in adata.obs["time"].unique():
    print(i)
    adata.obs["time"].replace([i],[count], inplace=True)
    count= count + 1
# adata.obs["time"].replace(['ResyncTime1', 'ResyncTime2', 'ResyncTime3', 'ResyncTime4', 'ResyncTime5',
#  'ResyncTime6'],[1, 2, 3, 4, 5, 6], inplace=True)
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
with plt.rc_context():
    sc.pl.embedding(
        adata,
        basis="X_draw_graph_fa",
        color=["time", annot_label],
        color_map="gnuplot",
        show=False)
    plt.tight_layout()
    plt.savefig(out_dir + "force_directed_graph.pdf", bbox_inches='tight')
    
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
print(gene_pool)

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
tp = tp.prepare(time_key="time")

# We solve one OT problem per time point pair, probabilistically matching 
# early to late cells.
print("Probalistically matching early to late cells")
tp = tp.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean")
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
    
# identfy ancestry of cells
print("identfy ancestry of cells")
# investigate which ancestry population a certain cell type has. 
# We do this by aggregating the transport matrix by cell type, using cell_transition()
# forward=True means we are plotting descendents, 
# forward=False means we are plotting ancestors
for i in adata.obs[annot_label].unique():
    print(i)
    clusters= adata.obs[annot_label].unique()
    for j in clusters:
        if i==j:
            pass
        else:
            key="transitions_"+i+"_"+j
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
                    fontsize=8, save=out_dir + key + ".png")

