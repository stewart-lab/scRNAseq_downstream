# Gamm_scRNAseq
Code for the analysis of the pig retinal organoid data from the David Gamm group.

Clone this repository:

```
git clone git@github.com:stewart-lab/Gamm_scRNAseq.git
```
You will also need anaconda or miniconda to install environments. To install miniconda:

https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html

## For normal single-cell processing without cross-species, visit our companion single-cell repo:

https://github.com/stewart-lab/scRNAseq_library

## Pre-processing for cross-species
R notebooks for processing the pig and human data. We take the orthologs across species, but process each species as a separate seurat object. NOTE: R notebooks can be run by opening up in R-studio or VS code, and running either by chunk, or running the whole thing at once.

First we need to install the single cell environment using conda and activate it
```
conda env create -f environment_scRNAseqbest.yml
conda activate scRNAseq_best
```
To run the following scripts- make sure to change your working directory, the GitHub directory with this repository, as well as the directory where your gene expression, cell, and gene matrices are in the .Rmd.

You will also need to update the directories for the orthologs from Ensemble and the meta data for the human reference in the .Rmd (both are provided in the data/ folder).

Variables:
```
WD <- working directory
GIT_DIR <- Github directory for Gamm_scRNAseq
ORTHOLOGS <- ortholog file
METADATA_REF <- meta data for human reference
```

Now run:
```
preprocess_crossspecies_Cowan.Rmd
preprocess_crossspecies_Reh.Rmd
```

## Seurat mapping
For mapping annotations across species, we used Seurat mapping.
For more details on Seurat mapping check out their webpage: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

First make sure the single cell environment we used previously is activated:
```
conda activate scRNAseq_best
```
Next, change the working directory, the GitHub directory with this repository, and directories where your pre-processed seurat objects are (for both human and pig) in the .Rmd file. Also update where the metadata files are for both human and pig (included in the data/ folder).

Variables:
```
WD <- working directory
GIT_DIR <- Github directory for Gamm_scRNAseq
METADATA_REF <- metadata for reference
METADATA_GAMM <- metadata for gamm data
```

Now run:
```
seurat_mapping_GAMM_Cowan.Rmd
seurat_mapping_GAMM_Reh.Rmd
```

## SC-comp
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

## Palantir and CellRank

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
