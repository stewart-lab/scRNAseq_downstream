# load packages

library(CellChat)
library(patchwork)
library(Seurat)
library(Matrix)
library(NMF)
library(ggalluvial)
library(reticulate)
options(stringsAsFactors = FALSE)

# load Seurat objects
# for seurat version < 5.0.0 normalized data: seurat_object[["RNA"]]@data # normalized data matrix
data_dir <- "/w5home/bmoore/Pierre_sc_zebrafish/"
setwd(data_dir)
# filename_list <- c("seurat_mouse_annot_4wkR.rds", 
#         "seurat_mouse_annot_6wkNR.rds", "seurat_mouse_annot_6wkR.rds")
filename_list <- read.table(file="filelist.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
filename_list <- filename_list$V1
print(filename_list)
# Create an empty list to store the imported objects
imported_objects <- list()
print("loading data")
# Loop through the filenames
for (i in seq_along(filename_list)) {
    # Extract the base filename without extension
    base_name <- tools::file_path_sans_ext(basename(filename_list[i]))
    # Import the file
    imported_object <- readRDS(filename_list[i])
    # Assign the imported object to the list with the base filename as the name
    imported_objects[[base_name]] <- imported_object
}
print(imported_objects)
object_names <- names(imported_objects)
seurat_vector <- unlist(imported_objects)

# get cell chat DB
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# subset DB
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# ID over-expression
id_overexpr <- function(cellchatobj){
    cellchatobj <- identifyOverExpressedGenes(cellchatobj)
    cellchatobj <- identifyOverExpressedInteractions(cellchatobj)
    return(cellchatobj)
}
# calculate cc communication prob
cc_prob <- function(cellchatobj, output, name){
        ptm = Sys.time()
        # options: # type set to trimean for fewer interactioins, also control for abundant cell type bias 
        cellchatobj <- computeCommunProb(cellchatobj, type = "triMean",population.size = TRUE)
        # filter few cells 
        cellchatobj <- filterCommunication(cellchatobj, min.cells = 1)
        # extract dfs
        df.net <- subsetCommunication(cellchatobj) # df of inferreed cc communications from ligand/receptors
        df.net.p <- subsetCommunication(cellchatobj, slot.name = "netP") # inferred at signal pathways
        # get signifcant pathways
        cellchatobj <- computeCommunProbPathway(cellchatobj)
        execution.time = Sys.time() - ptm
        print(paste0("communication porbability computation time: " ,as.numeric(execution.time, units = "secs")))
        # save
        print("saving dataframes")
        write.table(df.net,file=paste0(output, name, "_cellchat_df_net_ligand-recept.txt"), sep="\t", quote=FALSE)
        write.table(df.net.p,file=paste0(output, name, "_cellchat_df_net_signal_paths.txt"), sep="\t", quote=FALSE)
        return(cellchatobj)
}
# aggregate network
cc_aggregate <- function(cellchatobj, output, name){
    ptm = Sys.time()
    cellchatobj <- aggregateNet(cellchatobj)
    execution.time = Sys.time() - ptm
    print(paste0("aggregaing cell-cell network time: " ,as.numeric(execution.time, units = "secs")))
    # visualize as circle plot
    groupSize <- as.numeric(table(cellchatobj@idents))
    pdf(file=paste0(output, name, "_cc_interactions_celltype_circleplot.pdf"))
        par(mfrow = c(1,2), xpd=TRUE)
        netVisual_circle(cellchatobj@net$count, vertex.weight = groupSize, 
            weight.scale = T, label.edge= F, title.name = "Number of interactions")
        netVisual_circle(cellchatobj@net$weight, vertex.weight = groupSize, weight.scale = T, 
            label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    # save weight matrix
    write.table(cellchatobj@net$weight, file=paste0(output, name, "_cc_interactions_celltype_weight_matrix.txt"), 
        sep="\t", quote=FALSE, rownames=FALSE)
    mat <- cellchatobj@net$weight
    # visualize from each cell group
    pdf(file=paste0(output, name, "_cc_interactions_celltype_circleplot2.pdf"))
        par(mfrow = c(3,4), xpd=TRUE)
        for (i in 1:nrow(mat)) {
            mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
            mat2[i, ] <- mat[i, ]
            netVisual_circle(mat2, vertex.weight = groupSize, 
            weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
            }
    dev.off()
    return(cellchatobj)
}
# analyze signaling pathways
analyze_paths <- function(cellchatobj, vertex.receiver, output){
    # Access all the signaling pathways showing significant communications
    pathways.show.all <- cellchatobj@netP$pathways
    # make pathways folder
    path_dir <- paste0(output,"pathway_analysis")
    dir.create(path_dir, mode = "0777", showWarnings = FALSE)
    path_dir <- paste0(path_dir, "/")
    setwd(path_dir)
    # pathway analysis and figures
    for (i in 1:length(pathways.show.all)) {
            # Visualize communication network associated with both signaling pathway and individual L-R pairs
            print(pathways.show.all[i])
            netVisual(cellchatobj, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, 
                layout = "hierarchy",out.format ="pdf")
            # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
            gg <- netAnalysis_contribution(cellchatobj, signaling = pathways.show.all[i])
            ggsave(filename=paste0(path_dir, pathways.show.all[i], "_path_L-R_contribution.pdf"), 
                        plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
            # aggregate interactions
            pdf(file=paste0(path_dir, "cc_interactions_",pathways.show.all[i],"_circplot.pdf"))
            par(mfrow=c(1,1))
            netVisual_aggregate(cellchatobj, signaling = pathways.show.all[i], layout = "circle")
            dev.off()
                
            # Heatmap
            p2 <- netVisual_heatmap(cellchatobj, signaling = pathways.show.all[i], color.heatmap = "Reds")
            pdf(file=paste0(path_dir, "cc_interactions_",pathways.show.all[i],"_heatmap.pdf"))
            par(mfrow=c(1,1))
            print(p2)
            dev.off()
            # extract all sig L-R pairs and related sig genes for a pathway
            pairLR <- extractEnrichedLR(cellchatobj, signaling = pathways.show.all[i], geneLR.return = FALSE)
            LR.show <- pairLR[1,] # show one ligand-receptor pair
                
            # Circle plot
            pdf(file=paste0(path_dir, "top_L-Rpairs_",pathways.show.all[i],"_signaling.pdf"))
            netVisual_individual(cellchatobj, signaling = pathways.show.all[i], pairLR.use = LR.show, layout = "circle")
            dev.off()
            # bubble plot with sig interactions for a signaling pathway
            p4 <- netVisual_bubble(cellchatobj, sources.use = c(6), targets.use = c(4,12), 
                    signaling = pathways.show.all[i], remove.isolate = FALSE)
            pdf(file=paste0(path_dir,"sig_cc_L-Rpairs_",pathways.show.all[i],"_path_bubble_Tcells.pdf"), height = 11, width = 8.5)
            print(p4)
            dev.off()
            p5 <- netVisual_bubble(cellchatobj, sources.use = c(8), targets.use = c(6,5,3,7,13,14), 
                    signaling = pathways.show.all[i], remove.isolate = FALSE)
            pdf(file=paste0(path_dir,"sig_cc_L-Rpairs_",pathways.show.all[i],"_path_bubble_fibroblast.pdf"), height = 11, width = 8.5)
            print(p5)
            dev.off()  
            # plot gene expression for a given pathway
            p6 <- plotGeneExpression(cellchatobj, signaling = pathways.show.all[i], enriched.only = TRUE, type = "violin")
            pdf(file=paste0(path_dir,pathways.show.all[i],"_path_sig_geneExpr_violin.pdf"))
            print(p6)
            dev.off()            
        }
    setwd(output)
    return(cellchatobj)
}

# LR pair vizualization
LRpair_viz <- function(cellchatobj, output, name){
    # sig L-R pairs for interaction
    p1 <- netVisual_bubble(cellchat, sources.use = c(6), targets.use = c(4,12), remove.isolate = FALSE)
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_bubble_Tcells.pdf"), 
        height = 11, width = 8.5)
        print(p1)
    dev.off()
    # sig L-R pairs for interaction
    p2 <- netVisual_bubble(cellchat, sources.use = c(8), targets.use = c(6,5,3,7,13,14), remove.isolate = FALSE)
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_bubble_fibroblasts.pdf"), 
        height = 11, width = 8.5)
        print(p2)
    dev.off()
    # sig L-R pairs for interaction T cells
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_chord_Tcells.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchatobj, sources.use = c(6), targets.use = c(4,12), 
            lab.cex = 0.5,legend.pos.y = 30); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
            
    dev.off()
    
    # receiving
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_chord_Tcells_receive.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchatobj, sources.use = c(4,12), targets.use = c(6), 
            lab.cex = 0.5,legend.pos.y = 30); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
            
    dev.off()
    # sig L-R pairs for interaction fibroblasts
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_chord_fibroblasts.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchatobj, sources.use = c(8), targets.use = c(6,5,3,7,13,14), 
            legend.pos.x = 15); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    dev.off()
    # receiving
    pdf(file=paste0(output, name, "_sig_cc_L-Rpairs_celltypes_chord_fibroblasts_receive.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchatobj, sources.use = c(6,5,3,7,13,14), targets.use = c(8), 
            legend.pos.x = 15); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    dev.off()
    # show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
    pdf(file=paste0(output, name, "_sig_cc_pathways_chord_Tcells.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchat, sources.use = c(6), targets.use = c(4,12), 
            slot.name = "netP", legend.pos.x = 10); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    dev.off()
    pdf(file=paste0(output, name, "_sig_cc_pathways_chord_fibroblasts.pdf"))
        tryCatch({ result <- netVisual_chord_gene(cellchatobj, sources.use = c(8), targets.use = c(6,5,3,7,13,14), 
            slot.name = "netP", legend.pos.x = 10); print(result)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    dev.off()

    return(cellchatobj)
}

# net analysis
net_analysis <- function(cellchatobj, output, name){
    # Systems analysis of cell-cell communication network
    cellchatobj <- netAnalysis_computeCentrality(cellchatobj, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
    # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
    #Visualize dominant senders (sources) and receivers (targets)
    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    gg1 <- netAnalysis_signalingRole_scatter(cellchatobj)
    # save
    pdf(file=paste0(output, name, "_NetAnalysis_celltype_allpaths_2Dmap.pdf"))
        print(gg1)
    dev.off()
    # Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    ht1 <- netAnalysis_signalingRole_heatmap(cellchatobj, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchatobj, pattern = "incoming")
    pdf(file=paste0(output, name, "_NetAnalysis_allcelltype_allpaths_heatmap.pdf", 
        width = 10, height = 10))
        print(ht1 + ht2)
    dev.off()

    return(cellchatobj)
}

# Global communication patterns
global_cc <- function(cellchatobj, output, name, nPatternsOut, nPatternsIn){
    # set patterns for outgoing
    print("outgoing patterns")
    tryCatch({cellchatobj <- identifyCommunicationPatterns(cellchatobj, 
            pattern = "outgoing", k = nPatternsOut)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("outgoing pattern error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    # river plot
    print("plot")
    tryCatch({pr <- netAnalysis_river(cellchatobj, pattern = "outgoing")
            pdf(file=paste0(output, name, "_global_cc_patterns_outgoing.pdf"))
            print(pr)
            dev.off()},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    
    # dot plot
    tryCatch({pd <- netAnalysis_dot(cellchatobj, pattern = "outgoing")
            pdf(file=paste0(output, name, "_global_cc_patterns_outgoing_dotplot.pdf"))
            print(pd)
            dev.off()},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    
    # set patterns for incoming
    print("incoming patterns")
    tryCatch({cellchatobj <- identifyCommunicationPatterns(cellchatobj, 
                            pattern = "incoming", k = nPatternsIn)},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("incoming pattern error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    # river plot
    print("plot")
    tryCatch({pr2 <- netAnalysis_river(cellchatobj, pattern = "incoming")
            pdf(file=paste0(output, name, "global_cc_patterns_incoming.pdf"))
            print(pr2)
            dev.off()},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})
    
    # dot plot
    tryCatch({pd2 <- netAnalysis_dot(cellchatobj, pattern = "incoming")
            pdf(file=paste0(output, name, "global_cc_patterns_incoming_dotplot.pdf"))
            print(pd2)
            dev.off()},
            error= function(e) {an.error.occurred <<- TRUE
            message(paste("error? ",an.error.occurred))
            message("Here's the original error message:")
            message(conditionMessage(e))})

    return(cellchatobj)
}

# structure/ function analysis of signaling networks
funstr_analysis <- function(cellchatobj, output, name){
    # Identify signaling groups based on their functional similarity
    cellchatobj <- computeNetSimilarity(cellchatobj, type = "functional")
    #umap manifoldlearning
    cellchatobj <- netEmbedding(cellchatobj, type = "functional")
    # classification of signaling network
    cellchatobj <- netClustering(cellchatobj, type = "functional")
    # Visualization in 2D-space
    pf <- netVisual_embedding(cellchatobj, type = "functional", label.size = 3.5)
    pdf(file=paste0(output, name, "_pathway_functionality_2Dplot.pdf"))
    print(pf)
    dev.off()
    # Identify signaling groups based on structure similarity
    cellchatobj <- computeNetSimilarity(cellchatobj, type = "structural")
    cellchatobj <- netEmbedding(cellchatobj, type = "structural")
    cellchatobj <- netClustering(cellchatobj, type = "structural")
    #        Visualization in 2D-space
    ps <- netVisual_embedding(cellchatobj, type = "structural", label.size = 3.5)
    pdf(file=paste0(output, name, "_pathway_structural_2Dplot.pdf"))
    print(ps)
    dev.off()
    # zoom in
    pz <- netVisual_embeddingZoomIn(cellchatobj, type = "structural", nCol = 2)
    pdf(file=paste0(output, name, "_pathway_structural_2Dplot_zoomin.pdf"))
        print(pz)
    dev.off()

    return(cellchatobj)
}

# loop through to run cell chat
for(i in 1:length(object_names)){
        
        print(object_names[i])
        name <- object_names[i]

        #### make output dir  and read in data ####
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        output <- paste0(data_dir, name,"_cellchat_",timestamp)
        print(output)
        dir.create(output, mode = "0777", showWarnings = FALSE)
        output <- paste0(output, "/")
        seurat_obj <- seurat_vector[[name]]
        print(seurat_obj)

        #### make cell chat object ####
        cellchat <- createCellChat(object = seurat_obj, 
                    group.by = "CellType", assay = "RNA")
        print(cellchat)

        #### format cellchat obj ####
        # set idents
        cellchat <- setIdent(cellchat, ident.use = "CellType")
        levels(cellchat@idents)
        groupSize <- as.numeric(table(cellchat@idents))
        print(groupSize)
        # set the used database in the object
        cellchat@DB <- CellChatDB.use
        # preprocess expr data
        # subset the expression data of signaling genes for saving computation cost
        cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

        #### set maxSize for computation ####
        oopts <- options(future.globals.maxSize = 2.0 * 1e9) # 2GB
        on.exit(options(oopts))
        future::plan("multisession", workers = 4) # do parallel

        #### ID over-expression ####
        print("Identifying over-expressed genes and interactions")
        cellchat <- id_overexpr(cellchat)
        
        #### compute communication probability ####
        print("compute communication probability")
        cellchat <- cc_prob(cellchat, output, name)
        
        #### aggregate cell-cell network ####
        print("aggregaing cell-cell network")
        cellchat <- cc_aggregate(cellchat, output, name)
        
        #### pathway analysis ####
        print("doing individual pathway analysis")
        # check the order of cell identity to set suitable vertex.receiver
        levels(cellchat@idents)
        vertex.receiver = seq(1,4)
        # run pathway analysis
        cellchat <- analyze_paths(cellchat, vertex.receiver, output)
        
        #### visualize LR pairs ####
        print("visualize significant L-R pairs")
        cellchat <- LRpair_viz(cellchat, output, name)
        
        #### Network Analysis ####
        print("do Net analysis")
        cellchat <- net_analysis(cellchat, output, name)
        
        #### Global communication patterns ####
        print("select k for global communication patterns")
        # run selectK to infer the number of patterns.
        print("select k for outgoing")
        po <- selectK(cellchat, pattern = "outgoing")
        pdf(file=paste0(output, name, "_pattern_selection_outgoing.pdf"))
        print(po)
        dev.off()
        # incoming
        print("select k for ingoing")
        pi <- selectK(cellchat, pattern = "incoming")
        pdf(file=paste0(output, name, "_pattern_selection_incoming.pdf"))
        print(pi)
        dev.off()
        nPatternsIn = 5
        nPatternsOut = 4
        # where do Cophenetic and Silhouette values begin to drop suddenly?
        print("Running global communication patterns")
        cellchat <- global_cc(cellchat, output, name, nPatternsOut, nPatternsIn)
        
        #### functional and structural analysis ####
        print("manifold and classification learning analysis of signaling networks")
        # Manifold and classification learning analysis of signaling networks
        # quantify the similarity between all significant signaling pathways and 
        # then group them based on their cellular communication network similarity.
        # activate correct python env
        reticulate::use_python("/w5home/bmoore/.virtualenvs/cellchat/bin/python")
        # run analysis
        cellchat <- funstr_analysis(cellchat, output, name)

        #### SAVE ####
        print(paste0("saving object: ", name))
        saveRDS(cellchat, file=paste0(output, name, "_cellchat_obj.rds"))
        writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))

}
