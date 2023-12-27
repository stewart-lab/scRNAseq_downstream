# Gamm_scRNAseq
Code for the analysis of the pig retinal organoid data from David Gamm

Clone this repository:

```
git clone git@github.com:stewart-lab/Gamm_scRNAseq.git
```

## Pseudotime analysis
Cells are often in transition from one cell type to another, and pseudotime captures relationships between clusters, beginning with the least differentiated state to the most mature/terminal state(s).

## Palantir and CellRank

Palantir models trajectories of differentiated cells by treating cell fate as a probabilistic process and leverages entropy to measure cell plasticity along the trajectory. CellRank uses Palantir in it's pseudotime kernel, but can also use RNA velocity, similarity, cytotrace, time-series, or matebolic labeling to calculate trajectories. Here we use it with Palantir. Together they identify initial and terminal states, cell fate probabilities, and driver genes.

To run Palantir and CellRank you first have to have a Seurat object that is clustered and annotated (see above). 

Next convert your Seurat object to an anndata object (Note install the R requirements before running):

```
# before running, change the working directory and the input and output filenames in the script.

src/convert_seurat2anndata.R
```

- R requirements
```
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
library(SeuratDisk)
```

Once you have your anndata object, set up your python environment:

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
