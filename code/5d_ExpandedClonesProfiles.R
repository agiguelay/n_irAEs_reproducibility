##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#       STEP5d: Complementary work on the DE genes       #
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
library(ComplexHeatmap)
library(ggpubr)
library(AUCell)

# Parameters --------------------------------------------------------------------
dir.output <- "..." 
col.samples <- ...

# Data --------------------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI" &
                       FinalAnnotation == "CD4 CTL" & Top10_CD4_CTL != "No")
seurat.all <- subset(seurat.all, subset = FinalAnnotation == "CD4 CTL" & Top10_CD4_CTL != "No")
meta.good.cells <- seurat.all@meta.data

# Seurat workflow --------------------------------------------------------------

seurat.all <- FindVariableFeatures(seurat.all, selection.method = "vst", nfeatures = 2000)
genes.var <- VariableFeatures(seurat.all)
genes.to.remove <- grep("^TRAV|^TRBV|^TRBJ|^TRGV|^TRDV", genes.var, value = T)
genes.var <- setdiff(genes.var, genes.to.remove)

seurat.all <- ScaleData(seurat.all, verbose = FALSE, features = rownames(seurat.all))
seurat.all <- RunPCA(seurat.all, npcs = 30, verbose = FALSE, features = genes.var)
seurat.all <- RunUMAP(seurat.all, reduction = "pca", dims = 1:30, verbose = FALSE)

saveRDS(seurat.all, paste0(dir.output, "Objects/SeuratNonIntegrated_CD4CTL_Top10ExpandedCells.rds"))

pdf(paste0(dir.output, "Plots/CD4_CTLexp/UMAP_Seurat_Samples_CD4CTLexp.pdf"), width = 7, height = 5)
print(DimPlot(seurat.all, reduction = "umap", pt.size = 0.2, shuffle = T, group.by = "Sample", label = F,
                             label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
     scale_color_manual(values = col.samples) +
     theme(text = element_text(colour = "black", size = 15),
                     aspect.ratio = 1/1,
                     axis.line = element_line(color = "black", linewidth = 0.8),
                     axis.text = element_text(colour = "black", size = 15),
                     axis.ticks = element_line(color = "black", linewidth = 0.8)
               )
dev.off()

pdf(paste0(dir.output, "Plots/CD4_CTLexp/UMAP_Seurat_irAE_CD4CTLexp.pdf"), width = 7, height = 5)
print(DimPlot(seurat.all, reduction = "umap", pt.size = 0.2, shuffle = T, group.by = "irAE", label = F,
                              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
     scale_color_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
    theme(text = element_text(colour = "black", size = 15),
                     aspect.ratio = 1/1,
                    axis.line = element_line(color = "black", linewidth = 0.8),
                     axis.text = element_text(colour = "black", size = 15),
                    axis.ticks = element_line(color = "black", linewidth = 0.8)
               )
dev.off()


# PCA on differentially expressed genes ----------------------------------------
res <- read.table(paste0(dir.output, "Tables/DGE_Top10_CD4_CTL_ControlvsirAE_Post_PB.txt"), header = T)
res <- res[res$FDR <= 0.01 & abs(res$logFC) >= log2(1.5),]
genes <- rownames(res)
genes.to.remove <- grep("^TR", genes, value = T)
genes <- setdiff(genes, genes.to.remove)

seurat.all <- RunPCA(seurat.all, npcs = 30, verbose = FALSE, features = genes, reduction.name = "pca_sel") #new pca

pdf(paste0(dir.output, "Plots/CD4_CTLexp/PCA_Sel_Seurat_Samples_CD4CTLexp.pdf"), width = 7, height = 5)
print(DimPlot(seurat.all, reduction = "pca_sel", pt.size = 0.2, shuffle = T, group.by = "Sample", label = F,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = col.samples) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )

print(DimPlot(seurat.all, reduction = "pca_sel", pt.size = 0.2, shuffle = T, group.by = "irAE", label = F,
              label.size = 3, repel = TRUE, raster = F) + ggtitle("Query transferred labels")) +
  scale_color_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
dev.off()


## at the clonotype level
mat <- seurat.all@assays$RNA@scale.data
mat <- mat[genes,]
mat2 <- sapply(unique(seurat.all$clonotype_sample), 
               function(x) rowMeans(mat[, rownames(seurat.all@meta.data[seurat.all$clonotype_sample == x,])]))
pca <- prcomp(mat2)
df.pca <- data.frame(PCA1 = pca$rotation[,1], PCA2 = pca$rotation[,2], clonotype_sample = rownames(pca$rotation))
df.pca <- left_join(df.pca, meta.good.cells)

pdf(paste0(dir.output, "/Plots/CD4_CTLexp/PCA_Sel_Seurat_Samples_CD4CTLexp_Clono.pdf", width = 7, height = 5))
print(ggplot(df.pca, aes(x = PCA1, y = PCA2, color = irAE)) +
        geom_point() +
  scale_color_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
)
dev.off()

# Correlation between expansion and gene score ---------------------------------
genes.score <- list(irAEexp = intersect(genes, rownames(res[res$logFC < 0,])), 
                   CtrlExp = intersect(genes, rownames(res[res$logFC > 0,]))) # Calculate score
exprMatrix <- seurat.all@assays$RNA@counts
cells_AUC <- AUCell_run(exprMatrix, genes.score)
res_auc <- getAUC(cells_AUC)
res_auc <- as.data.frame(t(res_auc))

df <- cbind(meta.good.cells, res_auc)

pdf(paste0(dir.output, "Plots/CD4_CTLexp/Vln_irAEvsCtrl_CD4CTLexp_Scores_test.pdf"), width = 7, height = 5)
print(ggplot(df, aes(y = irAEexp, x = irAE, fill = irAE)) +
  geom_violin() +
  scale_fill_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
  stat_compare_means(method = "wilcox.test", paired = FALSE) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)))

print(ggplot(df, aes(y = CtrlExp, x = irAE, fill = irAE)) +
        geom_violin() +
        scale_fill_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)))

dev.off()

## Scores at the clonotype level
df.clono <- df %>% group_by(clonotype_sample) %>% 
  summarize(irAEexp = mean(irAEexp, na.rm = TRUE), CtrlExp = mean(CtrlExp, na.rm = TRUE)) %>%
  as.data.frame()
df.clono <- left_join(df.clono, meta.good.cells[!duplicated(meta.good.cells$clonotype_sample),])

pdf(paste0(dir.output, "Plots/CD4_CTLexp/Scatter_irAEvsCtrl_CD4CTLexp_Scores_Clono.pdf"), width = 7, height = 5)
print(ggplot(df.clono, aes(y = irAEexp, x = CtrlExp, col = irAE)) +
  geom_point() +
  scale_color_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
  theme(text = element_text(colour = "black", size = 15),
        aspect.ratio = 1/1,
        axis.line = element_line(color = "black", linewidth = 0.8),
        axis.text = element_text(colour = "black", size = 15),
        axis.ticks = element_line(color = "black", linewidth = 0.8)))

print(ggplot(df.clono, aes(y = irAEexp, x = CtrlExp, col = Sample)) +
        geom_point() +
        scale_color_manual(values =  col.samples) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)))
dev.off()


pdf(paste0(dir.output, "Plots/CD4_CTLexp/Vln_irAEvsCtrl_CD4CTLexp_Scores_Clono.pdf"), width = 7, height = 5)
print(ggplot(df.clono, aes(y = irAEexp, x = irAE, fill = irAE)) +
        geom_violin() +
        scale_fill_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)))

print(ggplot(df.clono, aes(y = CtrlExp, x = irAE, fill = irAE)) +
        geom_violin() +
        scale_fill_manual(values =  setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)))

dev.off()

pdf(paste0(dir.output, "Plots/CD4_CTLexp/BP_irAEvsCtrl_CD4CTLexp_Scores_Clono.pdf"), width = 7, height = 5)
print(ggplot(df.clono, aes(x = irAE, y = irAEexp, group.by = irAE)) +
        geom_boxplot(color = "black", linewidth = 1, outlier.shape = NA, fill = "white")  +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.2), dotsize = 1, stackratio = 0.5, aes(fill = Sample)) +
        scale_fill_manual(values = col.samples) +
        theme_classic() +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)
print(ggplot(df.clono, aes(x = irAE, y = CtrlExp, group.by = irAE)) +
        geom_boxplot(color = "black", linewidth = 1, outlier.shape = NA, fill = "white")  +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.2), dotsize = 1, stackratio = 0.5, aes(fill = Sample)) +
        scale_fill_manual(values = col.samples) +
        theme_classic() +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)
dev.off()

# UMAP Scores and genes --------------------------------------------------------
genes <- ...

for(g in genes){
  pdf(paste0(dir.output, "Plots/CD4_CTLexp/UMAP_", g, ".pdf"), width = 7, height = 5)
  print(FeaturePlot(seurat.all, reduction = "umap", pt.size = 0.2, features = g, raster = F) +
    ggtitle(g) +
    scale_color_gradientn(colors = c("#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
    theme(text = element_text(colour = "black", size = 15),
          aspect.ratio = 1/1,
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.text = element_text(colour = "black", size = 15),
          axis.ticks = element_line(color = "black", linewidth = 0.8)))
  dev.off()
}