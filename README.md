# scRNAseq_downstream

Code for the analysis of the pig retinal organoid data from the David Gamm group. Also contains other downstream applications to scRNAseq including different types of annotation using a reference (Seurat Mapping, scType), integration of datasets using Seurat, compositional analysis with scComp, pseudotime analysis, realtime analysis, and other processes.

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
For mapping annotations across species, we used Seurat mapping. Seurat mapping also uses a reference object/ dataset to predict annotations of clusters based on similarity of cell expression to the reference. It does this by doing a canonical correlation analysis to find "anchor" cells between the reference and query, then annotated clusters based on these anchor cells. 

For more details on Seurat mapping check out their webpage: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

First: **Make sure query and reference data were preprocessed the same way using the scRNAseq_library repo.**

Activate the single cell environment:
```
conda activate scRNAseq_new
```

Then update variables in the config file under the main header and seurat_mapping header.

Variables in config file:
```
"METHOD": "seurat_mapping
"seurat_mapping": {
    "DATA_DIR": "working_directory" # working directory where output goes,
    "REF.SEURAT": "ref.seurat.rds" # location (relative to DATA_DIR) and file name of reference seurat (from pre-processing),
    "QUERY.SEURAT": "query.seurat.rds" # location (relative to DATA_DIR) and file name of query seurat (from pre-processing),
    "get_metadata": {
      "get_meta": "TRUE" or "FALSE", # do you want to get metadata for object?
       "metadata_file1": "ref_metadata.txt" # location and filename of metadata for reference
       "metadata_file2": "query_metadata.txt" # location and filename of metadata for query
       "metadata_subset1": "d205" # if needed, what to subset the reference metadata by (ie sample name), else "NA"
       "metadata_subset2": "NA" # if needed, what to subset the query metadata by, else "NA"
    },
    "transfer_anchors": {
      "reduc.type": "cca", # type of reduction to use for anchors: "cca" or "pca"
      "query_manual_annot": "CellType_manual" # manual annotation label in query for comparisons
    },
    "visualize_and_subset_ref": {
      "groupby": "type", # label column in reference used for annotating
      "celltype_removal_list": ["AC1","AC2","T1/T3","T2","Midbrain","miG"] # list of cell types in reference to remove
    },
    "get_manual_comparison": { # order of prediction cell types ("rowvec") compared to manual cell types ("colvec")
    for proportion and cell count tables.
      "rowvec": ["BC","PR","iMG","Prog","Prog/Glia","AC","RGC","HC"],
      "colvec": ["Bipolar Cells","Rods","Cones","Rods - Muller Glia","Muller Glia","Muller Glia - Retinal Prog","Retinal Prog","Amacrine cells","unknown"]
    }
  }
```

Now run:
```
source run_downstream_toolkit.sh
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

To run, first activate the single cell environment:
```
conda activate scRNAseq_new
```

Then set variables under the "clustifyr" header. 

Variables in config file:
```
"METHOD":"clustifyr",
"clustifyr":{
    "DATA_DIR": "working_directory/", # working directory
    "REF.SEURAT": "../human_D205_subset_annot.rds", # location (relative to working dir) and file name of reference seurat (from pre-processing). If "NA", then only markers are used to annotate
    "QUERY.SEURAT": "../gamms2_cca_pred.rds", # location (relative to working dir) and file name of query seurat (from pre-processing)
    "cluster_name": "seurat_clusters", # cluster name to annotate
    "visualize_and_subset_ref": { # not used if no reference
      "groupby": "cell_type2", # label column in reference used for annotating
      "celltype_removal_list": ["AC-MC","CdBC-MC","L/M cone-rod","RPE"] # cell types not to be used for annotation
    },
    "score_and_plot_markers": {
      "known_markers": "True", # "True" if using known marker list, otherwise "False"
      "known_markers_path": "../../known_markers/Kims_retinal_markers.txt", # path where known markers are relative to working dir
      "cluster_type": "seurat_clusters", # which clusters to annotate
      "reduction": "umap" # dim reduction for visualization, ie. "umap", "pca"
    }
  }
```

Now run:
```
source run_downstream_toolkit.sh
```
Outputs:
* labeled Seurat object
* 2 query labeled annotation visualizations (one for marker list, one for reference)

## Annotation via marker lists

### GPT CellType
Automatic annotation of cell types using GPT-4 and differentially expressed marker genes. 

For more information see: https://winnie09.github.io/Wenpin_Hou/pages/gptcelltype.html

To run, first activate the single cell environment:
```
conda activate scRNAseq_new
```

Then set variables under the "celltypeGPT" header. 

Variables in config file:
```
"METHOD":"celltypeGPT"
"celltypeGPT":{
    "DATA_DIR": "working_directory/", # working directory,
    "seurat.obj": "seurat_obj_labeled.rds", # seurat object to annotate
    "openAI_key": "", # open AI key
    "cluster_name": "seurat_clusters2", # clusters to annotate
    "tissue_name": "retina" # if the tissue is known, otherwise NA
  }
```
Now run:
```
source run_downstream_toolkit.sh
```
Outputs:
* seurat object with annotations from GPT Celltype: seurat_celltypeGPT_annot.rds
* UMAP figure: celltypeGPT_umap.pdf

### scType
ScType a computational method for automated selection of marker genes based on scRNA-seq expression data. A cell type specificity score is assigned to each marker gene, and this is weighted by the gene expression matrix to ultimately get a scType score and cell type annotation for each cluster. See: https://github.com/IanevskiAleksandr/sc-type?tab=readme-ov-file

To run, first activate the single cell environment:
```
conda activate scRNAseq_new
```
Then set variables: 

Variables in config file (under sctype header):
```
"METHOD":"sctype"
"sctype":{
    "WD": "/w5home/bmoore/Brown_thymus/output_preprocess_20240910_095557/", # working dir
    "SEURAT_OBJ": "seurat_obj_labeled.rds", # clustered seurat object
    "db": "/w5home/bmoore/Brown_thymus/thymus_markers_YB_forScType.xlsx", # marker database (excel file with markers)
    "tissue": "thymus" # tissue type (should be present in marker list)
  }
```
Now run:
```
source run_downstream_toolkit.sh
```
Outputs:
* seurat_obj_labeled_sctype.rds --> seurat object with scType annotations
* sctype_umap.pdf --> umap of scType annotations
* sctype_bubbleplot.pdf --> bubble plot of all potential annotations per cluster

## Integration using Seurat
Integrate multiple Seurat objects. Objects are merged, then feature selection, scaling, and dimensionality reduction are performed. Next integration is done via canonical correlation analysis (CCA). After integration, clustering is performed and umap reduction is run again on the cca reduction, and integrated data are visualized. Finally layers are joined to then find differentially expressed genes across the integrated clusters.

To run, first activate the single cell environment:
```
conda activate scRNAseq_new
```

Variables in config file under "seurat_integration":
```
METHOD:"seurat_integration"
"seurat_integration":{
    "DATA_DIR":"working_dir/", # working directory
    "filename_list": ["seurat_obj_labeled_48h.rds","seurat_obj_labeled_72h.rds"], # list of seurat objects to be integrated with path relative to working directory
    "project_name": "GSE199571_merged_timepoints_with_fluidigm",
    "feature_selection": {
      "n_features": 2000, # number of genes to use
      "analysis_type": "Seurat" # must be "Seurat" for integrated data
    },
    "scale_data": {
      "vars.2.regress": "NA", # "NA" or "cell.cycle"
      "marker.path.s": "cell_cycle_vignette/cell_cycle_orthologs_s.genes.txt", # path to s cell cycle genes
      "marker.path.g2m": "cell_cycle_vignette/cell_cycle_orthologs_g2m.genes.txt" # path to g2m cell cycle genes
    },
    "run_and_visualize_pca": {
      "top_n_dims": 2, # number of pcs to visualize
      "heatmap_dims": 15, # number of pcs for heatmap
      "num_cells": 500, # number of cells to visualize on pca map
      "dims": 20, # number of dimensions to use for jack straw
      "num.replicate": "NA" # number of replicates to use for jack straw, or if "NA", jack straw not run
    },
    "run_umap": {
      "dims_umap": 20, # number of pcs to use for umap
      "umap.method": "umap-learn", # umap reduction method, "umap-learn" or "uwot"
      "umap.red": "pca" # reduction to use for umap, "pca" or "harmony"
    },
    "integration_type": "cca", # type of integration "cca", "rpca", "harmony", "mnn", or "scvi"
    "perform_clustering": { # clustering after integration
      "reduction": "integrated.cca", # use reduction according to integration type, ie. for "rpca" it is "integrated.rpca"
      "resolution": 0.5, # resolution for clusters
      "algorithm": "leiden", # louvain or leiden
      "dims_snn": 10, # dimensions of nearest neighbor graph to use
      "cluster_name": "cca_clusters" # cluster name 
    },
    "score_and_plot_markers": {
      "top_n_markers": 100, # Number of DE markers to consider
      "known_markers": "True", # if there are known markers- True or False
      "known_markers_path": "fluidigm_gup_expr_results/JackChu_markers.txt", # path with known markers file
      "cluster_type": "cca_clusters", # clusters to perform DE on
      "pairwise": "FALSE", # Compare each cluster combination pairwise? TRUE or FALSE
      "logFC_thresh": 0.25, # Cohen's D logFC threshold for a DE gene
      "auc_thresh": 0.5, # AUC cutoff (normally at least 0.5)
      "reduction": "umap.cca" # reduction to display umap, ie. shared space is umap.cca if cca was used
    },
    "process_known_markers":{
      "annot_type": "manual", # type of annotation
      "n_rank": 10 # lowest median logFC rank to consider for marker annotation
    }
  }
```
Now run.
```
source run_downstream_toolkit.sh
```
Outputs:
* combined seurat object
* pca and umap plots of integrated data
* DE genes for integrated clusters
* annotation if known marker list

## Cell type composition analysis

### SC-comp
To determine how the cell composition changed, we used sccomp, a Bayesian analysis that models changes in cell counts. For more information on sccomp: https://github.com/stemangiola/sccomp

First install and activate the sccomp environment:
```
conda env create -f environment_sccomp.yml
conda activate sccomp
```
Next update the .Rmd script and config with your working directory and the relative location of your **clustered** and **annotated** seurat objects.

Variables to set in config file under sccomp header:
```
"METHOD":"sccomp"
"sccomp":{
    "DATA_DIR": "/w5home/bmoore/scRNAseq/GAMM/human_data/reh_cellrep_2020/", # working dir
    "CROSS_SPECIES": "yes", # are you comparing across species (yes) or within species (no)
    "COMP_TYPE": "seuratmapping", # based on the annotation- did you use scpred, seuratmapping, or manual
    "annot_label1": "type", # annotation label for file 1
    "annot_label2": "predicted.id", # annotation label for file 2
    "idents1": "human_d205", # rename idents in file 1
    "idents2": "pig_day120", # rename idents in file 2
    "SEURAT.file1": "human_seurat_map.rds", # file 1 relative to working dir
    "SEURAT.file2": "query_seurat_pred.rds" # file 2 relative to working dir
  }
```
Now run.
```
source run_downstream_toolkit.sh
```
Outputs:
* combined seurat object
* proportion figure
* significance figures
* results tables

## Pseudotime analysis
Cells are often in transition from one cell type to another, and pseudotime captures relationships between clusters, beginning with the least differentiated state to the most mature/terminal state(s).

### Palantir and CellRank

Palantir models trajectories of differentiated cells by treating cell fate as a probabilistic process and leverages entropy to measure cell plasticity along the trajectory. CellRank uses Palantir in it's pseudotime kernel, but can also use RNA velocity, similarity, cytotrace, time-series, or matebolic labeling to calculate trajectories. Here we use it with Palantir. Together they identify initial and terminal states, cell fate probabilities, and driver genes.

To run Palantir and CellRank you first have to have a Seurat object that is clustered and annotated (see above). 

Next convert your Seurat object to an h5ad object that can be read into python. Modify working directory and object name in script
```
Rscript src/convert_seurat2anndata.R
```
Change to environments directory and install pseudotime environment
```
source install_pseudotime_env.sh
```
Activate a pseudotime environment
```
source pst_env/bin/activate
```

Modify config variables under pseudotime header:
```
"METHOD":"pseudotime"
"pseudotime":{
    "DATA_DIR": "/w5home/bmoore/scRNAseq/GAMM/", # working dir
    "ADATA_FILE": "seurat_obj_labeled.h5ad", # converted anndata file
    "NC": 10, # number of components that are used to find terminal cells. In general, lower for few terminal cell types, higher for many terminal cell types
    "METADATA": "manual_annot_metadata.txt", # if you need to add metadata
    "annot_label": "CellType1", # what annotation you want to analyze in pseudotime
    "gene_list": ["CNGA3","RCVRN", "CA10", "TUBB3","SOX2"], # genes that you want to visualize trajectories for
    "terminal_celltypes": ["Rod", "Cone", "Bipolar Cells", "Amacrine Cell"], # final cell types
    "start_celltype": "Retinal Prog", # progenitor or stem cell type
    "magic_impute": "FALSE", # impute gene expression data with MAGIC? TRUE or FALSE
    "trajectory":{ # pseudotime cell type trajectoires you want to visualize
      "cone": ["Retinal Prog", "unknown", "Pan PR", "Cone"],
      "rod": ["Retinal Prog", "unknown", "Pan PR", "Muller Glia-Rod", "Rod"],
      "MG": ["Retinal Prog", "unknown", "Muller Glia-Rod", "Muller Glia-Retinal Prog"],
      "bipolar": ["Retinal Prog", "unknown", "Bipolar Cells"],
      "amacrine": ["Retinal Prog", "unknown", "Amacrine Cell"]
    }
```
Now you are ready to run Palantir and CellRank. Change to src directory and run:
```
source run_downstream_toolkit.sh
```
Outputs:
* pseudotime_processed.h5ad object
* pseudotime figures including trajectories, palantir pseudotime, diffusion maps, transition maps, and terminal state probability figure
* gene figures including gene expression and gene trajectories over pseudotime

## Real time analysis
If you have actual time points in your data, you can do a real time analysis with optimal transport to see how the cells evolve over that time period.

### Realtime with Cell Rank and Moscot

First activate conda environment
```
conda activate realtime
```

Modify config variables under realtime header:
```
"METHOD":"realtime"
"realtime":{
    "DATA_DIR": "/w5home/bmoore/scRNAseq/LiFangChu/fluidigm_gup_expr_results/output_20231114_102348/", # working dir
    "ADATA_FILE": "clustered_seurat_obj.h5ad", # ann data object
    "time_label": "orig.ident", # label in object that designates time points
    "annot_label": "seurat_clusters" # annotation label to compare time points
  }
```
Run:
```
source run_downstream_toolkit.sh
```
Outputs:
* force_directed_graph.pdf --> force directed graph of time points and annotations
* prolif-apop_graph.pdf --> cell proliferation and apoptosis graph
* prior-post_growthrates.pdf --> prior and posterior cell growth rates
* cell_costs.pdf --> costs of source or target cell
* transitions_*.png --> figures identfiying descendents of celltypes

## Other processes

For other processes, use single cell environment:
```
conda activate scRNAseq_best
```
### Recluster
To re-cluster and re-annotate with updated cell types, use recluster-and-annotate.rmd:

Modify config variables under the recluster header:
```
"METHOD":"recluster"
"recluster":{
    "DATA_DIR": "/w5home/bmoore/Brown_thymus/output_preprocess_20240910_095557/PED165thymus/output_recluster_20240911_121040/",
    "SEURAT.FILE": "seurat_obj_labeled.rds",
    "perform_clustering": {
      "reduction": "pca",
      "resolution": 0.25,
      "algorithm": "leiden",
      "dims_snn": 10,
      "cluster_name": "seurat_clusters2"
    },
    "score_and_plot_markers": {
      "top_n_markers": 100,
      "known_markers": "True",
      "known_markers_path": "/w5home/bmoore/Brown_thymus/thymus_markers_3.txt",
      "cluster_type": "seurat_clusters2",
      "pairwise": "FALSE",
      "logFC_thresh": 0.25,
      "auc_thresh": 0.5,
      "reduction": "umap"
    },
    "process_known_markers":{
      "annot_type": "manual",
      "n_rank": 5
    }
```
Run:
```
source run_downstream_toolkit.sh
```
### Phate dimensionality reduction
To use the PHATE dimension reduction tool:
Modify config variables under "phate" header:
**options:**
* knn : Number of nearest neighbors (default: 5)
* decay : Alpha decay (default: 40)
* t : Number of times to power the operator (default: ‘auto’)
* gamma : Informational distance constant between -1 and 1 (default: 1)
```
"METHOD":"phate"
"phate":{
    "WD": "/w5home/bmoore/scRNAseq/LiFangChu/output_seuratintegrate_20240916_094725/",
    "SEURAT_OBJ": "seurat_obj_labeled.rds",
    "gene": "HES7",
    "knn":5, 
    "decay":50,
    "t":"auto",
    "gamma":1,
    "embed_key":"phate",
    "ANNOT": "CellType1"
  },
```
Run:
```
phate.Rmd
```
### Metadata
To output meta data to a text file, we used this script:
```
get_gamm_metadata.R
```
### Parse DE markers
To get the list of all DE genes and their annotation from each cluster for a seurat object (based on the output of a single cell analysis). Need a directory with the KnownDE.markers* files and Top100DEgenes_* files.
```
parse_markers.py <directory with DE gene output files>
```
### Make feature plots
Make feature plots for any list of genes to visualize expression in low dimensional space

Modify variables in config file under the header featureplots
```
"METHOD":"featureplots"
"featureplots":{
    "WD": "/w5home/bmoore/Brown_thymus/output_preprocess_20240910_095557/PED165thymus/output_recluster_20240911_121040/",
    "SEURAT_OBJ": "seurat_obj_labeled_phate.rds",
    "GENE_LIST": "../KM_genelist_thymus.txt",
    "ANNOT": "CellType1",
    "INPUT_NAME": "Thymus_0.5res_phate",
    "reduction": "phatek10d50"
  }
```
```
source run_downstream_toolkit.sh
```
### Subset cell cycle genes
Subset cell cycle genes by ortholog. Modify to input ortholog list and cell cycle gene list.
```
Rscript merge_orthologs_to_cellcycle_genes.R
```

# References: 

scPred paper: https://doi.org/10.1186/s13059-019-1862-5
**Note: scPRed is archived because it is not compatible with Seurat v5**

Seurat paper: https://doi.org/10.1016/j.cell.2019.05.031

sccomp paper: https://doi.org/10.1073/pnas.2203828120

Palantir paper: https://doi.org/10.1038/s41587-019-0068-4

CellRank2 paper: https://doi.org/10.1101/2023.07.19.549685

scType paper: https://doi.org/10.1038/s41467-022-28803-w

GPT CellType: https://doi.org/10.1038/s41592-024-02235-4

Clustifyr: https://doi.org/10.12688/f1000research.22969.2
