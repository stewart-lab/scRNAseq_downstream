
# sccomp analysis for cell type composition comparisons
# To run sccomp you need processed and labeled data sets that are annotated the same way

# Load libraries, set wd, and source functions
#### installations ####
# see scRNAseq_library/src/sccomp_env_setup.sh or
# use scRNAseq_library/environment_sccomp.yml to create an sccomp environment
## remember to set your Rterm in settings for remote to the conda path for sccomp!!!

### load libraries ###
print("load libraries and config")
library(reticulate)
#use_condaenv("sccomp", required=TRUE)
# devtools::install_github("Oshlack/speckle")
library(sccomp)
library(Seurat)
library(tidyverse)
library(loo)
library(reticulate)
library(devtools)
library(speckle)
library(ggplot2)
library(limma)
library(purrr)
library(jsonlite)
library(rmarkdown)
### set variables ###
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$sccomp$DATA_DIR
}
CROSS_SPECIES <- config$sccomp$CROSS_SPECIES # are you comparing across species (human to pig- then yes) or within species (no)
COMP_TYPE <- config$sccomp$COMP_TYPE # based on the annotation- did you use scpred, seuratmapping, or manual 
ANNOT_LABEL1 <- config$sccomp$annot_label1
ANNOT_LABEL2 <- config$sccomp$annot_label2
IDENTS1 <- config$sccomp$idents1
IDENTS2 <- config$sccomp$idents2
SEURAT.1 <- config$sccomp$SEURAT.file1
SEURAT.2 <- config$sccomp$SEURAT.file2

### set working directory and output ###
setwd(GIT_DIR)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_cell_composition_", timestamp)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(paste0(GIT_DIR, "/config.json"), file.path(output, "config.json"), 
        overwrite = TRUE)
print("sccomp version should be >= 1.4.0. package version: ")
packageVersion("sccomp") # should be >= 1.4.0

# Read in processed and labeled data
print("read in data")
seurat.obj.1 <- readRDS(file = paste0(DATA_DIR,SEURAT.1)) # reference
seurat.obj.2 <- readRDS(file = paste0(DATA_DIR,SEURAT.2)) # query

# check orig.ident and rename if necessary
print("check metadata")
colnames(seurat.obj.2@meta.data)
colnames(seurat.obj.1@meta.data)
table(seurat.obj.2$orig.ident)
table(seurat.obj.1$orig.ident)
# rename idents
Idents(object = seurat.obj.2) <- IDENTS2
Idents(object = seurat.obj.1) <- IDENTS1
# set as orig.ident
seurat.obj.1[["orig.ident"]] <- Idents(object = seurat.obj.1)
seurat.obj.2[["orig.ident"]] <- Idents(object = seurat.obj.2)
# check tables
table(seurat.obj.1$orig.ident)
table(seurat.obj.2$orig.ident)
# rename celltype column for reference
# reference cell type column
annot_ref<- seurat.obj.1[[ANNOT_LABEL1]]
seurat.obj.1$ref_type <- annot_ref
# check reference and query annotations
table(seurat.obj.1$ref_type)

# check location of annotation

if (CROSS_SPECIES=="yes"){
    #human
    table(annot_ref)
    #pig
    if (COMP_TYPE=="scpred"){
        table(seurat.obj.2$scpred_prediction)
    } else if (COMP_TYPE=="seuratmapping") {
       table(seurat.obj.2$predicted.id)
    } else {
        print("COMP_TYPE should be scpred or seurat mapping for cross-species")
    }
} else {
    if (COMP_TYPE=="scpred"){
        table(seurat.obj.1$scpred_prediction)
        table(seurat.obj.2$scpred_prediction)
    } else if (COMP_TYPE=="seuratmapping") {
       table(seurat.obj.1$predicted.id)
       table(seurat.obj.2$predicted.id)
    } else if (COMP_TYPE=="manual"){
        table(seurat.obj.1[[ANNOT_LABEL1]])
        table(seurat.obj.2[[ANNOT_LABEL2]])
    } else {
        print("COMP_TYPE should be scpred or seurat mapping or manual")
    }
}

# add identities to metadata

if (COMP_TYPE=="manual"){
    seurat.obj.1 <- AddMetaData(seurat.obj.1, metadata = seurat.obj.1[[ANNOT_LABEL1]], col.name= "CellType")
    table(seurat.obj.1$CellType)
    seurat.obj.2 <- AddMetaData(seurat.obj.2, metadata = seurat.obj.2[[ANNOT_LABEL2]], col.name= "CellType")
    table(seurat.obj.2$CellType)
}

# merge datasets
print("merge datasets by orig.ident")
seurat.combined <- merge(seurat.obj.1, y = seurat.obj.2, add.cell.ids = c("sample1", "sample2"), project = "combined")
head(colnames(seurat.combined))
table(seurat.combined$orig.ident)
unique(seurat.combined$orig.ident)
unique(sapply(X = strsplit(colnames(seurat.combined), split = "_"), FUN = "[", 1))
table(Idents(object = seurat.combined))

# set and rename identities in combined dataset
print("set identities in combined dataset")
# set identity

# stash idents
seurat.combined[["idents"]] <- Idents(object = seurat.combined)
unique(seurat.combined$idents)

# remove NAs
print("remove NAs")
# subset to remove NAs
seurat.combined <- subset(x= seurat.combined, idents != "NA")
seurat.combined <- subset(x= seurat.combined, orig.ident != "NA")

if (CROSS_SPECIES=="yes") {
    if (COMP_TYPE=="scpred") {
    # scpred
    seurat.pig <- subset(x= seurat.combined, scpred_prediction != "NA")
    unique(seurat.pig$scpred_prediction)
    } else if (COMP_TYPE=="seuratmapping") {
      # seurat mapping
    seurat.pig <- subset(x= seurat.combined, predicted.id != "NA")
    unique(seurat.pig$predicted.id)
    } else {
        print("COMP_TYPE should be scpred or seurat mapping for cross-species")
    }
    # human
    seurat.human <- subset(x= seurat.combined, ref_type != "NA")
    unique(seurat.human$ref_type)
}

# plot cell type proportions
print("plot proportions")
if (CROSS_SPECIES == "yes") {
    if (COMP_TYPE=="scpred") {
    # scpred
    pig.plot <- plotCellTypeProps(x = seurat.pig, clusters = seurat.pig$scpred_prediction, sample = seurat.pig$orig.ident)
    } else if (COMP_TYPE=="seuratmapping") {
    # seurat mapping
    pig.plot <- plotCellTypeProps(x = seurat.pig, clusters = seurat.pig$predicted.id, sample = seurat.pig$orig.ident)
    } else {
        print("COMP_TYPE should be scpred or seurat mapping for cross-species")
    }
    # human
    hum.plot <- plotCellTypeProps(x = seurat.human, clusters = seurat.human$ref_type, sample = seurat.human$orig.ident)
    # save plot
    pdf(paste0(output, "celltype_proportions.pdf"), width = 8, height = 6)
    print(pig.plot + hum.plot)
    dev.off()
} else {
    prop.plot <- plotCellTypeProps(x = seurat.combined, clusters = seurat.combined$CellType, sample = seurat.combined$idents)
    pdf(paste0(output, "celltype_proportions.pdf"), width = 8, height = 6)
    print(prop.plot)
    dev.off()
}

# combine cell types
print("combine cell types")
if (CROSS_SPECIES == "yes") {
    if (COMP_TYPE=="scpred") {
        New_idents <- c(seurat.human$ref_type,seurat.pig$scpred_prediction)
    } else if (COMP_TYPE=="seuratmapping") {
    New_idents <- c(seurat.human$ref_type,seurat.pig$predicted.id)
    } else {
        print("COMP_TYPE should be scpred or seurat mapping for cross-species")
    }
    seurat.combined@meta.data$"CellType" <- as.factor(New_idents)
}
table(seurat.combined$CellType)
# subset to remove NAs
seurat.combined <- subset(x= seurat.combined, CellType != "NA")

# save combined object
print("save combined object")
#saveRDS(seurat.combined, file = paste0(output,"seurat.combined.rds"))

# sccomp model composition
print("model composition")
# model composition
comp.tibble <- seurat.combined |>
  sccomp_estimate( 
    formula_composition = ~ idents, 
    .sample =  orig.ident, 
    .cell_group = CellType, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
# write table
comp.table <- as.data.frame(comp.tibble)
comp.table.1 <- comp.table[,1:7]
write.table(comp.table.1, file= paste0(output, "composition_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")
# get counts from a specific result:
celltype <- comp.table["12","CellType"]
pr.table <- as.data.frame(comp.table["12","count_data"])
celltype <- sub("/","_", celltype)

write.table(pr.table, file= paste0(output, celltype,"_composition_result_table.txt"), quote = F, col.names = TRUE, row.names= F, sep= "\t")


# compare models (Bayesian Anova)
print("compare models")
# need sccomp version 1.4 for this!
# Fit first model
model_with_factor_association <- 
  seurat.combined |>
  sccomp_estimate( 
    formula_composition = ~ idents, 
    .sample = orig.ident, 
    .cell_group = CellType, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
# Fit second model
model_without_association <- 
  seurat.combined |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    .sample = orig.ident, 
    .cell_group = CellType, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
# Compare models
comparison <- loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
comp_df <- as.data.frame(comparison)
write.table(comp_df, file= paste0(output, "model_comparison_table.txt"), 
            quote = F, col.names = TRUE, row.names= F, sep= "\t")
# is elpd_diff/se_diff > abs(5)?
c.res <- comp_df$elpd_diff[2]/comp_df$se_diff[2]
print(paste0("if ",c.res," is greater than abs(5), significant"))

# diffrential variability model
print("differential variability model")
res = 
  seurat.combined |>
  sccomp_estimate( 
    formula_composition = ~ idents,
    formula_variability = ~ idents,
    .sample = orig.ident,
    .cell_group = CellType,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
res.table <- as.data.frame(res)
res.table <- select(res.table, -c("count_data"))
write.table(res.table, file= paste0(output, "composition_x_var_result_table.txt"), 
            quote = F, col.names = TRUE, row.names= F, sep= "\t")

# test
test <- res |> sccomp_test()
test.table <- as.data.frame(test)
test.table <- select(test.table, -c("count_data"))
write.table(test.table, file= paste0(output, "composition_x_var_result_testtable.txt"), 
            quote = F, col.names = TRUE, row.names= F, sep= "\t")

# cell group statements

# res |> 
#    sccomp_proportional_fold_change(
#      formula_composition = ~  idents,
#      from =  IDENTS1,
#      to = IDENTS2
#     ) |> 
#   select(CellType, statement)



# plot
print("plot results")
# plot results
## for some reason boxplot not working
plots <- res |> sccomp_test() |> plot()
# group proportion boxplot 
# blue box is posterior predictive check and colors are 
# associated with significant associations for composition and/or variability
pdf(paste0(output, "signif_associations_boxplot.pdf"), width = 8, height = 6)
print(plots$boxplot)
# print(test |> sccomp_boxplot(factor="idents"))
dev.off()

# credible interval plot 1
# error bars represent 95% credible interval
pdf(paste0(output, "credible_intervals_1d.pdf"), width = 8, height = 6)
print(plots$credible_intervals_1D)
#print(test |> plot_1D_intervals())
dev.off()
# credible interval plot 2
pdf(paste0(output, "credible_intervals_2d.pdf"), width = 8, height = 6)
print(plots$credible_intervals_2D)
#print(test |> plot_2D_intervals())
dev.off()

# plot2 = model_with_factor_association |> sccomp_test() |> plot()
# # associated with significant associations for composition and/or variability
# pdf(paste0(output, "signif_associations_boxplot2.pdf"), width = 8, height = 6)
# print(plot2$boxplot)
# dev.off()
# # credible interval plot 1
# # error bars represent 95% credible interval
# pdf(paste0(output, "credible_intervals_1d2.pdf"), width = 8, height = 6)
# print(plot2$credible_intervals_1D)
# dev.off()
# # credible interval plot 2
# pdf(paste0(output, "credible_intervals_2d2.pdf"), width = 8, height = 6)
# print(plot2$credible_intervals_2D)
# dev.off()

# save session info
print("save session info")
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))
system(paste("chmod -R 777", output))