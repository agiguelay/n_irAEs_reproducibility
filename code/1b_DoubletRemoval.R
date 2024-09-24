##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#       STEP1b: Remove doublet with DoubletFinder        #
#                                                        #
#    On every individual Seurat object:                  #
#    Subsets the cells (QC)                              #
#    Runs standard Seurat workflow                       #
#    Runs DoubletFinder                                  #            
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################


rm(list = ls())

# Libraries --------------------------------------------------------------------

library(dplyr)
library(Seurat)
library(tidyverse)
library(SeuratData)
library(data.table)
library(DoubletFinder)
library(stringr)
library(ggplot2)

# Functions --------------------------------------------------------------------

# @param seurat.obj: Seurat object
# @param sample.ID: sample's name/ID
# @param dir.output: directory where all the outputs will be saved
# @param max.mtc: maximal value for mitochondrial content in a cell
# @param max.rbs: maximal value for ribosomal content in a cell
# @param min.features: minimal number of features detected in a cell
# @param min.umis: minimal number of UMIs detected in a cell
# @param max.umis: maximal number of UMIs detected in a cell (use a high number to avoid messing up with doublet removal)
# @param n.dim: number of PC to use
# Filter out low quality cells, runs standard Seurat's workflow and returns the
# seurat object

InitialDataProcessing <- function(seurat.obj, min.umis, max.umis, min.features, max.mtc, 
                                   max.rbs, col.palette, sample.ID,
                                   dir.output, n.dim = 30, ...){
  
  seurat.obj <- subset(seurat.obj, subset = nCount_RNA >= min.umis & 
                        nCount_RNA <= max.umis & nFeature_RNA >= min.features 
                        & percent.mt <= max.mtc & percent.rbs <= max.rbs)
  
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj, features = all.genes)
  
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:n.dim)
  seurat.obj <- FindClusters(seurat.obj, resolution = 0.5)
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:n.dim)
  
  df.plot <- data.frame(UMAP1 = seurat.obj@reductions$umap@cell.embeddings[,1],
                       UMAP2 = seurat.obj@reductions$umap@cell.embeddings[,2],
                       Cluster = seurat.obj$seurat_clusters)
  
  dir.create(paste0(dir.output, "Plots/DoubletRemoval/"), recursive = T)
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/UMAP_BeforeRemoval_", sample.ID, ".pdf"), width = 8, height = 5)
  print(ggplot(df.plot, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(alpha = 0.7, shape = 19, size = 0.5) + 
    scale_color_manual(values = sample(col.palette))+
    ggtitle(sample.ID) +
    theme_bw() +
    theme(text = element_text(colour = "black", size = 15), 
          aspect.ratio = 1/1,
          axis.line = element_line(color = "black", size = 0.8),
          axis.text = element_text(colour = "black", size = 15),
          axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  return(seurat.obj)
}



# from DoubletFinder but without plots
find.pK.cus = function(sweep.stats) {
  
  ## Implementation for data without ground-truth doublet classifications 
  '%ni%' <- Negate('%in%')
  if ("AUC" %ni% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=5))
    colnames(bc.mvn) <- c("ParamID","pK","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    # ## Plot for visual validation of BCmvn distribution
    # par(mar=rep(1,4))
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    # print(x)
    
    return(bc.mvn)
    
  }
  
  ## Implementation for data with ground-truth doublet classifications (e.g., MULTI-seq, CellHashing, Demuxlet, etc.)
  if ("AUC" %in% colnames(sweep.stats) == TRUE) {
    ## Initialize data structure for results storage
    bc.mvn <- as.data.frame(matrix(0L, nrow=length(unique(sweep.stats$pK)), ncol=6))
    colnames(bc.mvn) <- c("ParamID","pK","MeanAUC","MeanBC","VarBC","BCmetric")
    bc.mvn$pK <- unique(sweep.stats$pK)
    bc.mvn$ParamID <- 1:nrow(bc.mvn)
    
    ## Compute bimodality coefficient mean, variance, and BCmvn across pN-pK sweep results
    x <- 0
    for (i in unique(bc.mvn$pK)) {
      x <- x + 1
      ind <- which(sweep.stats$pK == i)
      bc.mvn$MeanAUC[x] <- mean(sweep.stats[ind, "AUC"])
      bc.mvn$MeanBC[x] <- mean(sweep.stats[ind, "BCreal"])
      bc.mvn$VarBC[x] <- sd(sweep.stats[ind, "BCreal"])^2
      bc.mvn$BCmetric[x] <- mean(sweep.stats[ind, "BCreal"])/(sd(sweep.stats[ind, "BCreal"])^2)
    }
    
    # ## Plot for visual validation of BCmvn distribution
    # par(mar=rep(1,4))
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, pch=18, col="black", cex=0.75,xlab=NA, ylab = NA)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$MeanAUC, col="black", lty=2)
    # par(new=TRUE)
    # x <- plot(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, pch=16, col="#41b6c4", cex=0.75)
    # axis(side=4)
    # x <- lines(x=bc.mvn$ParamID, y=bc.mvn$BCmetric, col="#41b6c4")
    # print(x)
    # 
    return(bc.mvn)
    
  }
}


# @param seurat.obj: Seurat object
# @param sample.ID: sample's name/ID
# @param dir.output: directory where all the outputs will be saved
# @param n.dim: number of PC to use
# @param dv.estim: estimated % of doublets (see 10x handbook)
# Filter out low quality cells, runs standard Seurat's workflow and returns the
# seurat object

DoubletRemoval <- function(seurat.obj, n.dim, db.estim, sample.ID, dir.output){ # Use DoubletFinder 
  
  ## pK identification
  sweep.res.list <- paramSweep_v3(seurat.obj, PCs = 1:n.dim, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK.cus(sweep.stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  pK <- bcmvn %>%
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  pK <- pK[[1]]
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/DoubletFinder_pKvalues_", sample.ID, '.pdf'), width = 6, height = 5)
  print(ggplot(bcmvn, aes(x = pK, y = BCmetric, group = 1)) +
    ggtitle(sample.ID) +
    geom_point(color = "deepskyblue4") + 
    geom_line(color = "deepskyblue4", size = 1) +
    geom_vline(xintercept = pK, linetype = "dashed", color = "orchid3", size = 1)+
    annotate(geom = "text", x = pK +0.04, y = max(bcmvn$BCmetric), label = paste0("pK = ", as.character(pK)), color = "orchid3", size = 6)+
    theme(text = element_text(colour = "black", size = 15), 
          axis.line = element_line(color = "black", size = 0.8),
          axis.text = element_text(colour = "black", size = 15),
          axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  ## Run DoubletFinder without taking into account homotypic doublet proportion
  nExp_poi <- round(db.estim*nrow(seurat.obj@meta.data))
  seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:n.dim, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
    ### UMAP projection
  col.db <- colnames(seurat.obj@meta.data)[grepl("DF.classification", colnames(seurat.obj@meta.data))]
  df.plot <- data.frame(UMAP1 = seurat.obj@reductions$umap@cell.embeddings[,1],
                       UMAP2 = seurat.obj@reductions$umap@cell.embeddings[,2],
                       DoubletStatus = seurat.obj@meta.data[, col.db])
  
  cat(paste("Without accounting for homotypic doublets: ", as.character(nExp_poi), " doublets \n" ))
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/UMAP_Doublets_", sample.ID, ".pdf"), width = 8, height = 5)
  print(ggplot(df.plot, aes(x = UMAP1, y = UMAP2, color = DoubletStatus)) +
          geom_point(alpha = 0.7, shape = 19, size = 0.5) + 
          scale_color_manual(values = setNames(c("violetred", "lightsteelblue"), c("Doublet", "Singlet")))+
          ggtitle(sample.ID) +
          theme_bw() +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1, 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
    ### QC plots
  df.plot$nCount_RNA <- seurat.obj$nCount_RNA
  df.plot$nFeature_RNA <- seurat.obj$nFeature_RNA
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Vln_nFeatures_Doublets_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nFeature_RNA, x = DoubletStatus)) +
          geom_violin(size = 0.8, fill = "#CC9933", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# detected features (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Vln_nUMIs_Doublets_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nCount_RNA, x = DoubletStatus)) +
          geom_violin(size = 0.8, fill = "#CCCC33", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# UMIs (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Scatter_nFeatures_vs_nUMIs_Doublets_", sample.ID, ".pdf"), width = 6, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nFeature_RNA, x = nCount_RNA, color = DoubletStatus)) +
          geom_point(size = 0.8, alpha = 0.5, shape  = 19) +
          scale_color_manual(values = setNames(c("violetred", "lightsteelblue"), c("Doublet", "Singlet")))+
          ylab("# detected features (RNA)") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1, 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  ## With Homotypic Doublet Proportion Estimate
  annotations <- seurat.obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  cat(paste("When accounting for homotypic doublets: ", as.character(nExp_poi.adj), " doublets \n" ))
  seurat.obj <- doubletFinder_v3(seurat.obj, PCs = 1:n.dim, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = str_replace(col.db, "DF.classifications", "pANN"), sct = FALSE)
  
  col.db2 <- colnames(seurat.obj@meta.data)[grepl("DF.classification", colnames(seurat.obj@meta.data))]
  col.db2 <- col.db2[grepl(as.character(nExp_poi.adj), col.db2)]
  df.plot$DoubletStatusHomotypic = seurat.obj@meta.data[, col.db2]
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/UMAP_DoubletsWithHomotypic_", sample.ID, ".pdf"), width = 8, height = 5)
  print(ggplot(df.plot, aes(x = UMAP1, y = UMAP2, color = DoubletStatusHomotypic)) +
          geom_point(alpha = 0.7, shape = 19, size = 0.5) + 
          scale_color_manual(values = setNames(c("violetred", "lightsteelblue"), c("Doublet", "Singlet")))+
          ggtitle(sample.ID) +
          theme_bw() +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1, 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Vln_nFeatures_DoubletsWithHomotypic_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nFeature_RNA, x = DoubletStatusHomotypic)) +
          geom_violin(size = 0.8, fill = "#CC9933", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# detected features (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Vln_nUMIs_DoubletsWithHomotypic_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nCount_RNA, x = DoubletStatusHomotypic)) +
          geom_violin(size = 0.8, fill = "#CCCC33", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# UMIs (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/DoubletRemoval/Scatter_nFeatures_vs_nUMIs_DoubletsWithHomotypic_", sample.ID, ".pdf"), width = 6, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = nFeature_RNA, x = nCount_RNA, color = DoubletStatusHomotypic)) +
          geom_point(size = 0.8, alpha = 0.5, shape  = 19) +
          scale_color_manual(values = setNames(c("violetred", "lightsteelblue"), c("Doublet", "Singlet")))+
          ylab("# detected features (RNA)") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  ## Gives a generic name
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = seurat.obj@meta.data[, col.db2], col.name = "DoubletStatus")
  return(seurat.obj)
}

