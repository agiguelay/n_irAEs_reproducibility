##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#     STEP3d: Pseudotime analyis (Slingshot) on CD4      #
#                                                        #
#    Integrates CD4 cells together and run Seurat        #
#    Runs Slingshot and find dynamic genes               #
#    Calculates the %irAE cells along pseudotime         #
#    Calculates the cell type composition along          #
#    pseudotime                                          #
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################

rm(list = ls())

# Libraries --------------------------------------------------------------------

library(Seurat)
library(slingshot)
library(ggplot2)
library(stringr)
library(SingleCellExperiment)
library(tradeSeq)
library(viridis)
library(ComplexHeatmap)
library(tradeSeq)
library(circlize)
library(RColorBrewer)
library(foreach)

# Parameters -------------------------------------------------------------------
dir.output <- "..."

col.palette <- c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                 "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                 "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")
col.cd4 <- setNames(col.palette[7:11], c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM"))

size.win <- 100 #half the length of the window for the sliding window analysis
n.perm <- 1000 #noumber of permutations for the sliding window analysis

# Load data and subsample ------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI" & FinalAnnotation %in% c("CD4 Naive","CD4 TCM", "CD4 TEM", "CD4 CTL")) # select the good cells, no proliferating cells because the clycling genes dominate too much

meta.good.cells <- seurat.all@meta.data

# Integration of CD4 T cells ---------------------------------------------------
files <- list.files(paste0(dir.output, "Objects/"), full.names = T)
files <- grep("SeuratObject", files, value = T)

meta.samples <- read.table("~/ag_ludwig/work/ICI_irAEs/SamplesMetadata_All.txt", sep = "\t", header = T)
meta.samples <- meta.samples[meta.samples$Dataset == "Maschmeyer" & meta.samples$Pre_Post == "Post-ICI",]
files <- grep(paste(meta.samples$Sample, collapse = "|"), files, value = T)

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
saveRDS(seurat.integrated, paste0(dir.output, "Objects/SeuratIntegrated_CD4_woTreg.rds"))

rm(our.anchors)
gc(reset = T)

DefaultAssay(seurat.integrated) <- "integrated"

# Seurat workflow --------------------------------------------------------------

seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE)
seurat.integrated <- RunPCA(seurat.integrated, npcs = 30, verbose = FALSE)
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)


seurat.integrated <- FindNeighbors(seurat.integrated)
seurat.integrated <- FindClusters(seurat.integrated, resolution = 0.3)
seurat.integrated <- RunUMAP(seurat.integrated, n.neighbors = 10, dims = 1:30, spread = 2, min.dist = 0.3)

seurat.integrated <- AddMetaData(seurat.integrated, meta.good.cells)
seurat.integrated <- AddMetaData(seurat.integrated, Idents(seurat.integrated), "ClusterInt")



meta.samples2 <- read.table("~/ag_ludwig/work/ICI_irAEs/SamplesMetadata_All.txt", sep = "\t", header = T, row.names = 1)
meta.samples2 <- meta.samples2[meta.samples2$Dataset == "Maschmeyer",]
meta.samples2 <- meta.samples2[order(meta.samples2$irAE),] #To have consistent colors between controls and irAE patients
col.samples <- setNames(col.palette[1:nrow(meta.samples2)], rownames(meta.samples2))

pdf(paste0(dir.output, "Plots/Pseudotime/UMAP_Seurat_Annotation_CD4_woTreg.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T, group.by = "FinalAnnotation", label = TRUE,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.cd4) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

pdf(paste0(dir.output, "Plots/Pseudotime/UMAP_Seurat_Clusters_CD4_woTreg.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T, group.by = "ClusterInt", label = TRUE,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.palette) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

pdf(paste0(dir.output, "Plots/Pseudotime/UMAP_Seurat_Samples_CD4_woTreg.pdf"), width = 12, height = 8)
print(DimPlot(seurat.integrated, reduction = "umap", pt.size = 0.2, shuffle = T, group.by = "Sample", label = F,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.samples) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()


# Run Slighshot ----------------------------------------------------------------
set.seed(1)
sce <- as.SingleCellExperiment(seurat.integrated)

sce <- slingshot(sce, clusterLabels = "ClusterInt", reducedDim = reducedDim(sce, 'PCA')[,1:30],
                allow.breaks = TRUE, start.clus = "0") # Starts in the naive cluster

dimred <- seurat.integrated@reductions$umap@cell.embeddings
clustering <- seurat.integrated$ClusterInt

pdf(paste0(dir.output, "Plots/Pseudotime/Lineage_PCA_CD4_woTreg.pdf"), width = 12, height = 8)
plot(reducedDims(sce)$PCA, col = col.palette[clustering], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black', type = "lineages")
dev.off()


colors <- viridis(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)] # Colors depending on pseudotime value

sds.new <- embedCurves(sce, dimred) # project the trajectories obtained with PCA on UMAP 

pdf(paste0(dir.output, "Plots/Pseudotime/PT_UMAP_CD4_woTreg.pdf"), width = 12, height = 8)
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
dev.off()

pdf(paste0(dir.output, "Plots/Pseudotime/PT_Linege_UMAP_CD4_woTreg.pdf"), width = 12, height = 8)
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sds.new), lwd = 3, type = "lineages")
dev.off()

seurat.integrated <- AddMetaData(seurat.integrated, setNames(sce$slingPseudotime_1, Cells(seurat.integrated)), "Pseudotime")
weights <- slingCurveWeights(sce)
pt <- slingPseudotime(sce)
pseudo.random <- sapply(1:nrow(pt), function(x) sample(pt[x,], 1, F, weights[x,]))
seurat.integrated <- AddMetaData(seurat.integrated, setNames(pseudo.random, Cells(seurat.integrated)), "PseudotimeRandomLin")
saveRDS(seurat.integrated, paste0(dir.output, "Objects/SeuratIntegrated_CD4_woTreg.rds"))

# Find genes correlating with pseudotime ---------------------------------------
counts <- seurat.integrated@assays$RNA@counts # add count assay for downstream analysis
counts(sce) <- counts

sce <- fitGAM(sce)
ATres <- associationTest(sce, contrastType="end", inverse="eigen") #Test for dynamic expression
write.table(ATres, paste0(dir.output, "Tables/DynamicGenes_CD4_woTreg.txt"), sep = "\t")

good.genes <- rownames(seurat.integrated)[apply(seurat.integrated@assays$RNA@counts, 1, function(x) sum(x > 0)) >= 50]
good.genes <- intersect(good.genes, rownames(seurat.integrated)[rowMeans(seurat.integrated@assays$RNA@data) >= 0.02])

topgenes <- ATres[order(ATres$pvalue), ]
topgenes <- topgenes[topgenes$pvalue <= 0.05,]
topgenes <- topgenes[intersect(rownames(topgenes), good.genes),]
topgenes <- rownames(topgenes)[1:250]
pst.ord <- order(seurat.integrated$Pseudotime, na.last = NA)
seurat.integrated <- ScaleData(seurat.integrated, verbose = FALSE, assay = "RNA")
mat <- seurat.integrated@assays$RNA@scale.data[topgenes, pst.ord]

ha <- HeatmapAnnotation(Pseudotime = seurat.integrated[,pst.ord]$Pseudotime, 
                       CellType = seurat.integrated[,pst.ord]$FinalAnnotation,
                       col = list(Pseudotime = colorRamp2(breaks = c(0, 20, 35, 45), colors = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")),
                                  CellType = col.cd4))

pdf(paste0(dir.output, "Plots/Pseudotime/HeatmapPT_CD4_woTreg.pdf"), width = 9, height = 17)

Heatmap(mat,
        show_column_names = F,
        column_order = colnames(mat),
        clustering_method_rows = "ward.D2",
        top_annotation = ha,
        row_names_gp = grid::gpar(fontsize = 5.5))
dev.off()

# Pseudotime distribution plots ------------------------------------------------

df.plot <- data.frame(Cell = Cells(seurat.integrated),
                     Condition = seurat.integrated$irAE,
                     Pseudotime = seurat.integrated$Pseudotime3,
                     CellType = seurat.integrated$FinalAnnotation,
                     Patient = seurat.integrated$Patient)

pdf(paste0(dir.output, "Plots/Pseudotime/DistributionPT3_CD4_woTreg.pdf"), width = 6, height = 4)
print(ggplot(df.plot, aes(x = Pseudotime, fill = Condition, color = Condition)) +
        geom_density(alpha= 0.5, aes(fill = Condition)) +
        theme_classic()+
        scale_color_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))
dev.off()

pdf(paste0(dir.output, "Plots/Pseudotime/DistributionPT3_CD4_woTreg_patient.pdf"), width = 6, height = 4)
print(ggplot(df.plot, aes(x = Pseudotime, color = Condition, groupby = Patient)) +
        geom_density(alpha= 0.5) +
        theme_classic()+
        scale_color_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))
dev.off()

# Sliding window PT analysis --------------------------------------------------
meta <- seurat.integrated@meta.data
meta <- meta[!is.na(meta$Pseudotime),]
meta <- meta[order(meta$Pseudotime),]
meta <- meta[, c("Pseudotime", "irAE")]



df <- data.frame(Pseudotime = meta[size.win:(nrow(meta)-size.win),]$PseudotimeRandomLin,
                NbiRAE = sapply(size.win:(nrow(meta)-size.win), function(x) sum(meta[(x - size.win):(x+size.win),]$irAE == "irAE")))
df$Prop_irAE <- (df$NbiRAE/(2*size.win+1))*100
df$Ratio <- df$NbiRAE/(2*size.win+1-df$NbiRAE)

pdf(paste0(dir.output, "Plots/Pseudotime/Prop_irAE_PT_CD4_woTreg.pdf"), width = 15, height = 5)
print(ggplot(df, aes(x = Pseudotime, y = Prop_irAE)) +
        geom_line(color = "#333399", linewidth = 1) +
        theme_classic() +
        geom_hline(yintercept = sum(meta$irAE == "irAE")/nrow(meta)*100, linetype = "dashed", color = "black") +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/2,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))

print(ggplot(df, aes(x = Pseudotime, y = Ratio)) +
        geom_line(color = "#333399", linewidth = 1) +
        theme_classic() +
        geom_hline(yintercept = sum(meta$irAE == "irAE")/sum(meta$irAE == "Control"), linetype = "dashed", color = "black") +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/2,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))

dev.off()

# Add some permutation tests --------------------------------------------------

prop.perm = foreach(i = 1:n.perm, .combine = "cbind", .inorder = F) %do% {
  meta.perm <- meta
  meta.perm$irAE <- sample(meta.perm$irAE)
  NbiRAE <- sapply(size.win:(nrow(meta.perm)-size.win), function(x) sum(meta.perm[(x - size.win):(x+size.win),]$irAE == "irAE"))
  Prop_irAE <- (NbiRAE/(2*size.win+1))*100
  out <- data.frame(Proportion = Prop_irAE, row.names = 1:nrow(df))
  out
}

write.table(prop.perm, paste0(dir.output, "Tables/PermutationPT_CD4.txt"))

df$LowQuantile <- apply(prop.perm, 1, function(x) quantile(x, 0.025))
df$HighQuantile <- apply(prop.perm, 1, function(x) quantile(x, 0.975))

pdf(paste0(dir.output, "Plots/Pseudotime/PropiRAPTRd_CD4_woTreg_Permutations_5.pdf"), width = 15, height = 5)
print(ggplot(df, aes(x = Pseudotime, y = Prop_irAE)) +
        geom_ribbon(aes(x = Pseudotime, ymax = HighQuantile, ymin = LowQuantile), fill = "#9999CC", alpha = .5) +
        geom_line(aes(y = HighQuantile), colour = '#666699') +
        geom_line(aes(y = LowQuantile), colour = '#666699')+
        theme_classic() +
        geom_hline(yintercept = sum(meta$irAE == "irAE")/nrow(meta)*100, linetype = "dashed", color = "black") +
        geom_line(color = "#333366", linewidth = 1) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/2,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))
dev.off()



## Cell Type proportion ----
meta <- seurat.integrated@meta.data
meta <- meta[!is.na(meta$Pseudotime),]
meta <- meta[order(meta$Pseudotime),]
meta <- meta[, c("Pseudotime", "FinalAnnotation")]

df2 <- data.frame(Pseudotime = meta[size.win:(nrow(meta)-size.win),]$Pseudotime,
                NbCTL = sapply(size.win:(nrow(meta)-size.win), function(x) sum(meta[(x - size.win):(x+size.win),]$FinalAnnotation == "CD4 CTL")),
                NbTEM = sapply(size.win:(nrow(meta)-size.win), function(x) sum(meta[(x - size.win):(x+size.win),]$FinalAnnotation == "CD4 TEM")),
                NbTCM = sapply(size.win:(nrow(meta)-size.win), function(x) sum(meta[(x - size.win):(x+size.win),]$FinalAnnotation == "CD4 TCM")),
                NbNaive = sapply(size.win:(nrow(meta)-size.win), function(x) sum(meta[(x - size.win):(x+size.win),]$FinalAnnotation == "CD4 Naive")))

df2$`CD4 CTL` <- (df2$NbCTL/(2*size.win+1))*100
df2$`CD4 TEM` <- (df2$NbTEM/(2*size.win+1))*100
df2$`CD4 TCM` <- (df2$NbTCM/(2*size.win+1))*100
df2$`CD4 Naive` <- (df2$NbNaive/(2*size.win+1))*100

df3 <- stack(df2[, 6:9])
colnames(df3) <- c("Proportion", "CellType")
df3$Pseudotime <- rep(df2$Pseudotime, 4)

pdf(paste0(dir.output, "Plots/Pseudotime/Prop_Celltype_PTRd_CD4_woTreg.pdf"), width = 15, height = 5)
print(ggplot(df3, aes(x = Pseudotime, y = Proportion, color = CellType)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = col.cd4) +
        theme_classic() +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/2,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))

dev.off()