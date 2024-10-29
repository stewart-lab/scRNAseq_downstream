## scRNAseq env setup

conda env create -n scRNAseq_new --file scRNAseq_new.yml

## in R:
# install.packages("remotes")
# remotes::install_github('satijalab/seurat-wrappers')
# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github('satijalab/azimuth', ref = 'master')
# install.packages("openxlsx", dependencies = TRUE)
# install.packages("openai")
# remotes::install_github("Winnie09/GPTCelltype")
