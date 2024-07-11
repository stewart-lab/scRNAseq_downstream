# sccomp environment

conda create -n sccomp r-essentials r-base
conda activate sccomp
conda install -c conda-forge r-reticulate
conda install -c bioconda r-seurat
conda install -c bioconda bioconductor-sccomp
conda install -c bioconda bioconductor-limma
conda install -c conda-forge r-devtools
conda install -c bioconda bioconductor-edger
conda install -c bioconda bioconductor-biocparallel

# in R do: devtools::install_github("Oshlack/speckle")