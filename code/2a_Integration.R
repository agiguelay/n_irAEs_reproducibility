##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#       STEP2a: Integrate the individual samples         #
#                                                        #
#    On every individual Seurat object:                  #
#    Removes doublets, calculates a cell cycle score     #
#    and find variable features                          #            
#    Integrates with Seurat's method                     #
#    Runs Seurat's workflow on the integrated object     #
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------
library(Seurat)
library(foreach)
library(stringr)
library(ggplot2)
library(data.table)


# Parameters -------------------------------------------------------------------
dir.output <- "..."

col.palette <- c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")


# Load seurat object, subset the cells, calculate cell cycle score and ---------
# identify variable genes
files <- list.files(paste0(dir.output, "Objects/"), full.names = T)
files <- grep("SeuratObject", files, value = T)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

libs <- NULL
list.obj <- list()
for(f in files){
  lib <- str_match(f, "SeuratObject_\\s*(.*?)\\s*.rds")[2]
  print(lib)
  libs <- c(libs, lib)
  seurat.obj <- readRDS(f)
  
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- subset(seurat.obj, subset = DoubletStatus == "Singlet")
  print(dim(seurat.obj))
  
  seurat.obj$Sample <- rep(lib, ncol(seurat.obj))
  seurat.obj <- CellCycleScoring(seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  seurat.obj <- RenameCells(seurat.obj, new.names = paste0(lib, "_", Cells(seurat.obj)))
  
  seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)
  seurat.obj <- FindVariableFeatures(seurat.obj, verbose = FALSE)
  saveRDS(seurat.obj, file <- paste0(dir.output, "Objects/SeuratObject_", lib, ".rds"))
  list.obj[[lib]] <- seurat.obj

}


# Seurat's integration ---------------------------------------------------------
features <- SelectIntegrationFeatures(object.list = list.obj)

list.obj <- lapply(list.obj, function(x) {
  x = subset(x, features = features)
  x = ScaleData(x, features = features, verbose = FALSE) #On RNA assay
  x = RunPCA(x, features = features, verbose = FALSE)
})

our.anchors <- FindIntegrationAnchors(object.list = list.obj,
                                     dims = 1:50, anchor.features = features,
                                     reduction = "rpca")

rm(list.obj)
gc(reset = T)

seurat.integrated <- IntegrateData(anchorset = our.anchors, dims = 1:50)

rm(our.anchors)
gc(reset = T)

DefaultAssay(seurat.integrated) <- "integrated"


# Seurat's standard workflow ---------------------------------------------------
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca",
                             dims = 1:30, verbose = FALSE)

dir.create(paste0(dir.output, "Plots/Integration/"), recursive = T)

pdf(paste0(dir.output, "Plots/Integration/UMAP_SeuratIntegration_Sample.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", group.by = "Sample", 
              pt.size = 0.7, shuffle = T, raster = F) +
        scale_color_manual(values = col.palette) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))
dev.off()

seurat.integrated <- FindNeighbors(seurat.integrated, dims = 1:30)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.8)

pdf(paste0(dir.output, "Plots/Integration/UMAP_SeuratIntegration_Clusters.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", group.by = "seurat_clusters",
              pt.size = 0.7, shuffle = T, raster = F) +
        scale_color_manual(values = col.palette) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))
dev.off()

clusters <- data.frame(Cluster = seurat.integrated$seurat_clusters, Cell = rownames(seurat.integrated@meta.data), Sample = seurat.integrated$Sample,
                      UMAP1 = Embeddings(seurat.integrated, reduction = "umap")[, 1], UMAP2 = Embeddings(seurat.integrated, reduction = "umap")[, 2])
write.table(clusters, paste0(dir.output, "Tables/Clustering_integration.txt"), sep = "\t")