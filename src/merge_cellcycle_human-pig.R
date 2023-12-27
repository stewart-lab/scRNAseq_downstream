# get orthologs and cell cycle genes
orthologs<- read.table("human_ref/Human_Pig_Biomart_Filtered_mod.txt", header=TRUE, sep="\t")
cell.cycle<- read.table("regev_lab_cell_cycle_genes.txt", header=TRUE, sep="\t")
# replace colname so its the same as the ortholog name
colnames(cell.cycle)[1] <- "human.gene.name"
# merge 
cell.cycle.orthos <- merge(cell.cycle,orthologs,by="human.gene.name") 
# write out
write.table(cell.cycle.orthos, file="cell_cycle_orthologs.txt",col.names=TRUE, sep="\t",quote=F )
# get s/g2m phases from seurat
library(Seurat)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
s.genes<- as.data.frame(s.genes)
colnames(s.genes)[1] <- "human.gene.name"
cell.cycle.orthos.s.genes <- merge(s.genes,cell.cycle.orthos,by="human.gene.name") 

g2m.genes <- cc.genes$g2m.genes
g2m.genes <- as.data.frame(g2m.genes)
colnames(g2m.genes)[1] <- "human.gene.name"
cell.cycle.orthos.g2m.genes <- merge(g2m.genes,cell.cycle.orthos,by="human.gene.name")

# write out s and g2m genes
write.table(cell.cycle.orthos.s.genes, file="cell_cycle_orthologs_s.genes.txt",col.names=TRUE, 
            sep="\t",quote=F )
write.table(cell.cycle.orthos.g2m.genes, file="cell_cycle_orthologs_g2m.genes.txt",col.names=TRUE, 
            sep="\t",quote=F )
