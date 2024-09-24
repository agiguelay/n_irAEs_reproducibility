##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#        STEP4b: Performs UMAP on Tcell compartment      #
#                                                        #
#    Integrates T cells together and run Seurat          #
#    Projects on a UMAP clonal informations              #          
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(stringr)
library(SingleCellExperiment)
library(tradeSeq)
library(viridis)
library(ComplexHeatmap)
library(tradeSeq)
library(foreach)

# Parameters -------------------------------------------------------------------
dir.output = "..."
col.ct = ... # For the cell types
col.samples = ... # For the samples

# Load data and subset ---------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
t.cells <- c("gdT", "CD4 TEM", "CD4 TCM","CD8 TEM","CD4 CTL", "CD8 Naive", 
            "CD8 TCM", "dnT", "Treg", "CD4 Proliferating", "MAIT", "CD8 Proliferating", 
            "CD4 Naive")
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI" & FinalAnnotation %in% t.cells) # select the good cells
meta.good.cells <- seurat.all@meta.data

# Integration of T cells -------------------------------------------------------
files <- list.files(paste0(dir.output, "Objects", full.names = T))
files <- grep("SeuratObject", files, value = T)

libs <- NULL
list.obj <- list()
for(f in files){
  lib <- str_match(f, "SeuratObject_\\s*(.*?)\\s*.rds")[2]
  print(lib)
  libs <- c(libs, lib)
  seurat.obj <- readRDS(f)
  
  DefaultAssay(seurat.obj) <- "RNA"
  cells.to.keep <- intersect(Cells(seurat.obj), rownames(meta.good.cells))
  seurat.obj <- subset(seurat.obj, cells = cells.to.keep)
  print(dim(seurat.obj))
  
  list.obj[[lib]] <- seurat.obj
}

features <- SelectIntegrationFeatures(object.list = list.obj)

list.obj <- lapply(list.obj, function(x) {
  
  x <- subset(x, features = features)
  x <- ScaleData(x, features = features, verbose = FALSE) #On RNA assay
  x <- RunPCA(x, features = features, verbose = FALSE)
})

our.anchors <- FindIntegrationAnchors(object.list = list.obj,
                                     dims = 1:50, anchor.features = features, reduction = "rpca")

rm(list.obj)
gc(reset = T)

seurat.integrated <- IntegrateData(anchorset = our.anchors, dims = 1:50)

rm(our.anchors)
gc(reset = T)

DefaultAssay(seurat.integrated) <- "integrated"

# Seurat workflow --------------------------------------------------------------

seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)


seurat.integrated <- FindNeighbors(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, n.neighbors = 10, dims = 1:30, spread = 2, min.dist = 0.3)

seurat.integrated <- AddMetaData(seurat.integrated, meta.good.cells)

# UMAP projection --------------------------------------------------------------

pdf(paste0(dir.output, "Plots/VDJ/UMAP_Seurat_Annotation_Tcells.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T,
              group.by = "FinalAnnotation", label = FALSE,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.ct) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()


pdf(paste0(dir.output, "Plots/VDJ/UMAP_Seurat_TopCloneCD4_Tcells.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T,
              group.by = "Top10_CD4", label = F,
              label.size = 3, repel = TRUE, raster = F)) +
  scale_color_manual(values = c(col.samples, setNames("gray85", "No"))) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

pdf(paste0(dir.output, "Plots/VDJ/UMAP_Seurat_TopCloneCD8_Tcells.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T,
              group.by = "Top10_CD8", label = F,
              label.size = 3, repel = TRUE, raster = F)) +
  scale_color_manual(values = c(col.samples, setNames("gray85", "No"))) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()
