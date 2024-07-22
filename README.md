# scRNAseq_downstream

Code for the analysis of the pig retinal organoid data from the David Gamm group. Also contains other downstream applications to scRNAseq including different types of annotation using a reference (Seurat Mapping, scPred), integration of datasets using Seurat, compositional analysis with scComp, and pseudotime analysis.

Clone this repository:

```
git clone git@github.com:stewart-lab/scRNAseq_downstream.git
```
You will also need anaconda or miniconda to install environments. To install miniconda:

https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html

## For single-cell pre-processing, visit our companion single-cell repo:

https://github.com/stewart-lab/scRNAseq_library

## Annotation via a reference

### Seurat mapping
For mapping annotations across species, we used Seurat mapping. Seurat mapping also uses a reference object/ dataset to predict annotations of clusters based on similairty of cell expression to the reference. It does this by doing a canonical correlation analysis to find "anchor" cells between the reference and query, then annotated clusters based on these anchor cells. 

For more details on Seurat mapping check out their webpage: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

First: **Make sure query and reference data were preprocessed the same way using the scRNAseq_library repo.**

Next activate the single cell environment we used previously is activated:
```
conda activate scRNAseq_best
```
Next, change the working directory, the GitHub directory with this repository, and directories where your pre-processed seurat objects are (for both human and pig) in the .Rmd file. Also update where the metadata files are for both human and pig (included in the data/ folder) in the config file.

Variables in .Rmd file:
```
WD <- working directory
GIT_DIR <- Github directory for scRNAseq_downstream
REF.SEURAT <- location and file name of reference seurat (from pre-processing)
QUERY.SEURAT <- location and file name of query seurat (from pre-processing)
```
Variables in config file:
```
in "get_metadata":
  "metadata_file1": location of metadata for reference
  "metadata_file2": location of metadata for query
  "metadata_subset1": if needed, what to subset the reference metadata by (ie sample name), else "NA"
  "metadata_subset2": if needed, what to subset the query metadata by, else "NA"
in "transfer_anchors":
    "reduc.type": type of reduction to use for anchors: "cca" or "pca"
    "query_manual_annot": manual annotation label in query for comparisons, i.e. "CellType_manual"
in "visualize_and_subset_ref":
    "groupby": label column in reference used for annotating, i.e. "type"
    "celltype_removal_list": list of cell types in reference to remove, i.e. ["AC2","miG","T2"]
in "get_manual_comparison": order of prediction cell types ("rowvec") compared to manual cell types ("colvec")
    for proportion and cell count tables.
    "rowvec": ["CdBC","ChBC","RBC","rod","L/M cone","MC","AC","HC"],
    "colvec": ["Bipolar Cells","Rods","Cones","Rods - Muller Glia","Muller Glia","Muller Glia - Retinal Prog",
"Retinal Prog","Amacrine cells","unknown"]
```

Now run:
```
Rscript -e "rmarkdown::render('seurat_mapping.Rmd')"
```

Outputs:
* annotated query Seurat object
* Subsetted reference object
* Umaps of transferred labels to query in query space and reference space
* Umap feature plots showing prediction score of each cell type
* Proportion and cell count tables comparing prediction to manual annotation
* Heatmap of proportion comparison to manual annotation

### ClustifyR
ClustifyR uses either a marker list of genes or a reference (or both!) to annotate a query object. For a marker list, Clustifyr annotates based on percent cells expressed in a given cluster for a given marker. It can also use enrichment tests. For a reference, Clustifyr performs a Spearman's correlation between the reference expression matrix and query expression matrix to find maximally correlated clusters, annotating with the highest correlated reference cluster (above a threshold).

For more information on ClustifyR see: https://www.bioconductor.org/packages/release/bioc/vignettes/clustifyr/inst/doc/clustifyr.html

To run, first set variables:

Variables in .Rmd file:
```
WD <- working directory
GIT_DIR <- Github directory for scRNAseq_downstream
REF.SEURAT <- location and file name of reference seurat (from pre-processing)
QUERY.SEURAT <- location and file name of query seurat (from pre-processing)
```
Variables in config file:
```
 in "score_and_plot_markers":
    "known_markers": "True" if using known marker list, otherwise "False"
    "known_markers_path": path where known markers are relative to working dir
    "cluster_type": which clusters to annotate
    "reduction": dim reduction for visulaization, ie. "umap", "pca"
in visualize_and_subset_ref:
    "groupby" : label column in reference used for annotating, i.e. "cell_type2"
```

Now run:
```
Rscript -e "rmarkdown::render('clustifyr.Rmd')"
```
Outputs:
* labeled Seurat object
* 2 query labeled annotation visualizations (one for marker list, one for reference)

## Integration using Seurat
Integrate multiple Seurat objects. Objects are merged, then feature selection, scaling, and dimensionality reduction are performed. Next integration is done via canonical correlation analysis (CCA). After integration, clustering is performed and umap reduction is run again on the cca reduction, and integrated data are visualized. Finally layers are joined to then find differentially expressed genes across the integrated clusters.

Variables in .Rmd file:
```
WD <- working directory
GIT_DIR <- Github directory for Gamm_scRNAseq
filename_list <- list of seurat objects to be integrated with path relative to working directory. i.e. c("output/object1.rds","output/object2.rds")
```

Variables in config file:
```
in "feature_selection":
    "n_features": number of genes to use, i.e. 2000
    "analysis_type": must be "Seurat" for integrated data
in "scale_data": 
    "vars.2.regress": "NA" or "cell.cycle"
    "marker.path.s": path to s cell cycle genes, ie. "../cell_cycle_vignette/cell_cycle_orthologs_s.genes.txt",
    "marker.path.g2m": path to g2m cell cycle genes, ie. "../cell_cycle_vignette/cell_cycle_orthologs_g2m.genes.txt"
in "run_and_visualize_pca": 
    "top_n_dims": number of pcs to visualize, i.e. 2
    "heatmap_dims": number of pcs for heatmap, i.e. 15
    "num_cells": number of cellls to visualize on pca map, i.e. 500
    "dims": number of dimensions to use for jack straw, i.e. 20
    "num.replicate": number of replicates to use for jack straw, i.e. 100. if "NA", jack straw not run
in "run_umap": 
    "dims_umap": number of pcs to use for umap, i.e. 20
    "umap.method": umap reduction method, "umap-learn" or "uwot"
    "umap.red": reduction to use for umap, "pca" or "harmony"
in "perform_clustering": (performed after integration)
    "reduction": "integrated.cca"
    "resolution": 0.5,
    "algorithm": "leiden"
    "dims_snn": 10
    "cluster_name": "cca_clusters"
in "score_and_plot_markers": 
    "top_n_markers": number of DE markers to consider, i.e. 100
    "known_markers": if known markers, "True" or "False"
    "known_markers_path": path to known markers
    "cluster_type": cluster type to annotate, i.e. "cca_clusters"
    "pairwise": perform pairwise DE between clusters, "TRUE" or "FALSE"
    "logFC_thresh": 0.25
    "auc_thresh": 0.5
"process_known_markers":
   "annot_type": "manual"
   "n_rank": lowest median rank to consider for marker annotation, i.e. 10
```

## Cell type composition analysis

### SC-comp
To determine how the cell composition changed, we used sccomp, a Bayesian analysis that models changes in cell counts. For more information on sccomp: https://github.com/stemangiola/sccomp

First install and activate the sccomp environment:
```
conda env create -f environment_sccomp.yml
conda activate sccomp
```
Next update the .Rmd script with your working directory and the relative location of your **clustered** and **annotated** seurat objects.
Variables to set in .Rmd:
```
BETHS_GAMM_DATA_DIR <- "/w5home/bmoore/scRNAseq/GAMM/" # main directory
CROSS_SPECIES <- "no" # are you comparing across species (human to pig- then yes) or within species (no)?
COMP_TYPE <- "manual" # based on the annotation- did you use scpred, seuratmapping, or manual
```
Now run.
```
sccomp.Rmd
```

## Pseudotime analysis
Cells are often in transition from one cell type to another, and pseudotime captures relationships between clusters, beginning with the least differentiated state to the most mature/terminal state(s).

### Palantir and CellRank

Palantir models trajectories of differentiated cells by treating cell fate as a probabilistic process and leverages entropy to measure cell plasticity along the trajectory. CellRank uses Palantir in it's pseudotime kernel, but can also use RNA velocity, similarity, cytotrace, time-series, or matebolic labeling to calculate trajectories. Here we use it with Palantir. Together they identify initial and terminal states, cell fate probabilities, and driver genes.

To run Palantir and CellRank you first have to have a Seurat object that is clustered and annotated (see above). 

Next convert your Seurat object to an h5ad object that can be read into python. 
**Note** install the R requirements before running:

- R requirements
```
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
library(SeuratDisk)
```

**Before running**, change the working directory and the input and output filenames in the script.

Variables:
```
WD <- working directory
SEURAT_OBJ <- seurat object name
METADATA_GAMM <- gamm metadata
OUTPUT_name <- output file name for the h5ad object that will be created
```
Now run.
```
Rscript src/convert_seurat2anndata.R
```

Once you have your h5ad object, set up your python environment:

```
# set up with pseudotime_requirements.txt

source src/install_pseudotime_env.sh

# activate

source pst_venv/bin/activate
```

Now you are ready to run Palantir and CellRank

- open src/pseudotime_GAMM.ipynb in vscode and update the following variables:

```
# variables
DATA_DIR = "your_directory" # directory where anndata object is
ADATA_FILE = "gamms2_cca_pred.h5ad" # the name of your h5ad (anndata) file
ANNOT_TYPE = "manual" # the type of annotation used, options are: "seurat_map", "clustifyr", "manual"
CROSS_SPECIES = "TRUE" # is this a cross-species annotation? "TRUE" or "FALSE"
NC = 8 # number of components that are used to find terminal cells. In general, lower for few terminal cell types, higher for many terminal cell types
```

- now run src/pseudotime_GAMM2.ipynb
- figures will be saved to the data directory

## Other processes

To modify the initial annotation with updated cell types, we used this R-markdown script:
```
reannotate_manual.rmd
```
To output meta data to a text file, we used this script:
```
get_gamm_metadata.R
```
To get the list of all DE genes and their annotation from each cluster for a seurat object (based on the output of a single cell analysis):
```
parse_markers.py
```
To get histograms of prediction scores:
```
predscore_stats.rmd
```
# References: 

scPred paper: https://doi.org/10.1186/s13059-019-1862-5

Seurat paper: https://doi.org/10.1016/j.cell.2019.05.031

sccomp paper: https://doi.org/10.1073/pnas.2203828120

Palantir paper: https://doi.org/10.1038/s41587-019-0068-4

CellRank2 paper: https://doi.org/10.1101/2023.07.19.549685
