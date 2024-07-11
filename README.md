# scRNAseq_downstream

Code for the analysis of the pig retinal organoid data from the David Gamm group. Also contains other downstream applications to scRNAseq including different types of annotation using a reference (Seurat Mapping, scPred), integration of datasets using Seurat, compositional analysis with scComp, and pseudotime analysis.

Clone this repository:

```
git clone git@github.com:stewart-lab/Gamm_scRNAseq.git
```
You will also need anaconda or miniconda to install environments. To install miniconda:

https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html

## For single-cell pre-processing with or without cross-species, visit our companion single-cell repo:

https://github.com/stewart-lab/scRNAseq_library

## Annotation via a reference

### Seurat mapping
For mapping annotations across species, we used Seurat mapping.Seurat mapping also uses a reference object/ dataset to predict annotations of clusters based on similairty of cell expression to the reference. It does this by doing a canonical correlation analysis to find "anchor" cells between the reference and query, then annotated clusters based on these anchor cells.

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
GIT_DIR <- Github directory for Gamm_scRNAseq
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
```

Now run:
```
Rscript -e "rmarkdown::render('seurat_mapping.Rmd')"
```
### scPred

scPred uses a reference object/ dataset to predict annotations of clusters in query data based on similarity of cell expression to the reference. The default algorithm is SVM-radial, however many different models/algorithms can be applied from the caret package: https://topepo.github.io/caret/available-models.html. 

The introduction vignette for scPred can be found here: https://powellgenomicslab.github.io/scPred/articles/introduction.html

To run, your query and reference data must first be processed the same way. 

Next run the scPred script with your preprocessed query and reference objects. scPred will divide the reference into training and testing objects, where the model is trained on the training set, and then applied to the testing set. Watch for cell types that don't predict well in the test set- this may mean the model for that cell type isn't good, and you can try a different one. After a final model is built (you can have different algorithms for each cell type if you want), then you can apply to the query data. To run, update with your reference and query objects, you can also specify different algorithms/models after first running your training data with the SVM-radial algorithm.

```
src/scPred_GAMM.Rmd
```

- R requirements
```
library(devtools)
devtools::install_github("powellgenomicslab/scPred")
library(scPred)
library(Seurat)
library(magrittr)
library(harmony)
library(rmarkdown)
library(jsonlite)
library(purrr)
library(scran)
library(patchwork)
library(dplyr)
library(reticulate)
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
