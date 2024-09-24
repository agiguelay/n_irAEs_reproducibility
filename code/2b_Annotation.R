##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#       STEP2b: Cell type annotation with Azimuth        #
#                                                        #
#    On every individual Seurat object:                  #
#    Runs Azimuth                                        #
#    Merge seurat objects and add the umap + clustering  #
#    calculated on the integrated object (2a)            #          
#    Check some markers/QC/prediction score              #
#    Readjust the annotation                             #
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
library(Azimuth)

# Parameters -------------------------------------------------------------------
dir.output <- "..."

col.palette <- c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                 "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                 "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")


# Run Azimuth ------------------------------------------------------------------
files <- list.files(paste0(dir.output,"Objects/"), full.names = T)
files <- grep("SeuratObject", files, value = T)

libs <- NULL
list.obj <- list()
for(f in files){
  lib <- str_match(f, "SeuratObject_\\s*(.*?)\\s*.rds")[2]
  print(lib)
  libs <- c(libs, lib)
  seurat.obj <- readRDS(f)
  
  DefaultAssay(seurat.obj) <- "RNA"
  
  # Transfer label with azimuth 
  seurat.obj <- RunAzimuth(seurat.obj, reference = "pbmcref", assay = "RNA")
  saveRDS(seurat.obj, file = paste0(dir.output, "Objects/SeuratObject_", lib, ".rds"))

  # seurat.obj[["prediction.score.celltype.l1"]] <- NULL
  # seurat.obj[["prediction.score.celltype.l2"]] <- NULL
  # seurat.obj[["prediction.score.celltype.l3"]] <- NULL
  # seurat.obj[["refAssay"]] <- NULL

  list.obj[[lib]] <- seurat.obj
  
}

# Merge all libraries ----------------------------------------------------------
seurat.all <- merge(x = list.obj[[1]], y = list.obj[2:22])
seurat.all@meta.data <- seurat.all@meta.data[,1:26]

# Add embeddings and clustering performed on the integrated object -------------
umap <- read.table(paste0(dir.output, "Tables/Clustering_integration.txt"), sep = "\t", header = T)

seurat.all$SeuratClusterIntegration <- umap[Cells(seurat.all),]$Cluster
seurat.all[["FinalUMAP"]] <- CreateDimReducObject(embeddings = data.matrix(umap[Cells(seurat.all), c("UMAP1", "UMAP2")]), key = "UMAP_", assay = DefaultAssay(seurat.all))
saveRDS(seurat.all, paste0(dir.output, "Objects/SeuratFinal_All.rds"))

col.palette = c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")

col.celltypes <- setNames(col.palette[1:length(unique(seurat.all$predicted.celltype.l2))], sort(unique(seurat.all$predicted.celltype.l2)))

pdf(paste0(dir.output, "Plots/Integration/UMAP_Seurat_AzimuthAnnotation.pdf"), width = 12, height = 8)
print(DimPlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, shuffle = T, group.by = "predicted.celltype.l2", label = TRUE,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.celltypes) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

# Check mapping and prediction scores ------------------------------------------

pdf(paste0(dir.output, "Plots/Integration/UMAP_Seurat_PredictionScore.pdf"), width = 8, height = 8)
print(FeaturePlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, features = c("predicted.celltype.l2.score"), raster = F) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1, 
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)

print(FeaturePlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, features = c("mapping.score"), raster = F) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)
dev.off()

# Check Quality control metrics ------------------------------------------------

pdf(paste0(dir.output, "Plots/Integration/UMAP_Seurat_Integration_QC.pdf"), width = 8, height = 8)
print(FeaturePlot(seurat.all, features = c("percent.mt"), reduction = "FinalUMAP", pt.size = 0.2, raster = T) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)

print(FeaturePlot(seurat.all, features = c("percent.rbs"), reduction = "FinalUMAP", pt.size = 0.2, raster = F) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)

print(FeaturePlot(seurat.all, features = c("nCount_RNA"), reduction = "FinalUMAP", pt.size = 0.2, raster = F) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)

print(FeaturePlot(seurat.all, features = c("nFeature_RNA"), reduction = "FinalUMAP", pt.size = 0.2, raster = F) +
        scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)
dev.off()


# Readjust azimuth annotation --------------------------------------------------

seurat.all$FinalAnnotation = seurat.all$predicted.celltype.l2
seurat.all$FinalAnnotation[seurat.all$SeuratClusterIntegration %in% c(21, 27)] = "Low QC cells"
seurat.all$FinalAnnotation[seurat.all$SeuratClusterIntegration == 14] = "Treg"
seurat.all$FinalAnnotation[seurat.all$SeuratClusterIntegration == 10] = "gdT"
seurat.all$FinalAnnotation[seurat.all$SeuratClusterIntegration == 17] = "CD4 CTL"

pdf(paste0(dir.output, "Plots/Integration/UMAP_Seurat_FinalAnnotation.pdf"), width = 12, height = 8)
print(DimPlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, shuffle = T, group.by = "FinalAnnotation", label = TRUE,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
        scale_color_manual(values = c(col.celltypes, setNames("gray85", "Low QC cells"))) +
        theme(text = element_text(colour = "black", size = 15),
             aspect.ratio = 1/1,
             axis.line = element_line(color = "black", linewidth = 0.8),
             axis.text = element_text(colour = "black", size = 15),
             axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

table(seurat.all$FinalAnnotation, seurat.all$predicted.celltype.l2) #check
saveRDS(seurat.all, "~/ag_ludwig/work/ICI_irAEs/Results/Objects/SeuratFinal_All.rds")

png("~/ag_ludwig/work/ICI_irAEs/Results/Plots/Integration/UMAP_Seurat_FinalAnnotation.png", width = 12, height = 8, units = "in", res = 600)
print(DimPlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, shuffle = T, group.by = "FinalAnnotation", label = F, raster = F)+
        scale_alpha_manual(values = rep(1, length(unique(seurat.all$FinalAnnotation)))) +
        scale_color_manual(values = c(col.celltypes, setNames("gray85", "Low QC cells"))) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        )
)
dev.off()

# Generates a dotplot with markers ---------------------------------------------
Idents(seurat.all) <- "FinalAnnotation"

ct <- c("CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg",
       "CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM",
       "MAIT", "gdT", "dnT", "ILC",
       "B naive", "B intermediate", "B memory", "Plasmablast",
       "NK", "NK_CD56bright", "NK Proliferating",
       "CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "ASDC", "pDC", 
       "Eryth", "Platelet", "HSPC", "Low QC cells")

Idents(seurat.all) <- factor(seurat.all$FinalAnnotation, levels = rev(ct))
markers <- c("CD3D", "CD4", "CD8A", "IL2RA", "FOXP3", "IL7R", "CD27", "GZMB", 
             "GZMA", "GZMK", "CCR7", "CD28","SELL", "TCF1", "TCF7", "MS4A1", 
             "CD79A", "NCAM1", "CD14", "FCGR3A", "LYZ", "CD68", "CLEC4C", "MKI67")

pdf(paste0(dir.output, "Plots/Integration/DotPlot_Seurat_Integration_Markers_1.pdf"), width = 10, height = 6)
DotPlot(
  object = seurat.all,
  assay = "RNA",
  features = markers,
  cols = "PuOr",
  col.min = -1.5,
  col.max = 1.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()