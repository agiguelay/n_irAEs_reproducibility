##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#           STEP5b: Heatmap of differentially            #
#         expressed genes (clonotype pseudobulk)         #
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


# Parameters --------------------------------------------------------------------
dir.output <- "..." 
col.samples <- ...


# Load data --------------------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI")

# For CD4 CTL ------------------------------------------------------------------

# Select good genes
res <- read.table(paste0(dir.output, "Tables/DGE_Top10_CD4_CTL_ControlvsirAE_Post_PB.txt"), header = T)
res <- res[res$FDR <= 0.01 & abs(res$logFC) >= log2(1.5),]
small.seurat <- subset(seurat.all, subset = FinalAnnotation == "CD4 CTL" & Top10_CD4_CTL != "No")
small.seurat <- ScaleData(small.seurat, features = rownames(small.seurat))
mat <- small.seurat@assays$RNA@scale.data
mat <- mat[rownames(res),]

genes.to.remove <- grep("^TR", rownames(mat), value = T)
mat <- mat[setdiff(rownames(mat), genes.to.remove),]

order.patient <- small.seurat@meta.data[order(small.seurat$irAE),]
order.patient <- unique(order.patient$Sample)
small.seurat$Sample <- factor(small.seurat$Sample, levels = order.patient)

# Pseudobulk: average per clonotype

mat2 <- sapply(unique(small.seurat$clonotype_sample), function(x) 
  rowMeans(mat[, rownames(small.seurat@meta.data[small.seurat$clonotype_sample == x,])]))
sample <- sapply(colnames(mat2), function(x) strsplit(x, "clonotype([0-9]+)_")[[1]][2])

# Heatmap 
ha <- HeatmapAnnotation(Condition = meta.samples[sample,]$irAE,
                        Sample = sample,
                        col = list(Condition = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")),
                                   Sample = col.samples))

sample <- sample[order(match(sample, order.patient))]

pdf(paste0(dir.output, "Plots/DGE/TopClono/HeatmapAverage_Post_Top10_CD4_CTL_PB.pdf"), height = 5, width = 6)
Heatmap(data.matrix(mat2), show_column_names = F, top_annotation = ha,
        clustering_method_rows = "ward.D2",
        row_names_gp = gpar(fontsize = 6),
        column_order = names(sample))
dev.off()

# Select genes for OR (5c)

mat.sel <- small.seurat@assays$RNA@data
mat.sel <- mat.sel[setdiff(rownames(res), genes.to.remove),]
mat.sel <- sapply(unique(small.seurat$clonotype_sample), 
                 function(x) rowMeans(mat.sel[, rownames(small.seurat@meta.data[small.seurat$clonotype_sample == x,])]))

tab.samples <- table(sample)
samples.solo <- names(tab.samples[tab.samples == 1])
clono.solo <- names(sample[sample %in% samples.solo])
samples.multi <- setdiff(unique(sample), samples.solo)

mat.sel2 <- cbind(sapply(samples.multi, function(x) # patient agregation
  rowSums(mat.sel[, names(sample[sample == x])])), mat.sel[, clono.solo])
sel.genes <- intersect(rownames(mat.sel2[apply(mat.sel2, 1, function(x) sum(x > 0)) >= 3,]), #in at least 3 patients
                      rownames(mat.sel[apply(mat.sel, 1, function(x) sum(x >=0.2)) >= 5,])) #minimum value in at least 5 clonotypes (i.e. 5%)
length(sel.genes)

# For CD8 TEM ------------------------------------------------------------------

# Select good genes
res <- read.table(paste0(dir.output, "Tables/DGE_Top10_CD8_TEM_ControlvsirAE_Post_PB.txt"), header = T)
res <- res[res$FDR <= 0.01 & abs(res$logFC) >= log2(1.5),]
small.seurat <- subset(seurat.all, subset = FinalAnnotation == "CD8 TEM" & Top10_CD8_TEM != "No")
small.seurat <- ScaleData(small.seurat, features = rownames(small.seurat))
mat <- small.seurat@assays$RNA@scale.data
mat <- mat[rownames(res),]

genes.to.remove <- grep("^TR", rownames(mat), value = T)
mat <- mat[setdiff(rownames(mat), genes.to.remove),]

order.patient <- small.seurat@meta.data[order(small.seurat$irAE),]
order.patient <- unique(order.patient$Sample)
small.seurat$Sample <- factor(small.seurat$Sample, levels = order.patient)

# Pseudobulk: average per clonotype
mat2 <- sapply(unique(small.seurat$clonotype_sample), function(x) rowMeans(mat[, rownames(small.seurat@meta.data[small.seurat$clonotype_sample == x,])]))

sample <- sapply(colnames(mat2), function(x) strsplit(x, "clonotype([0-9]+)_")[[1]][2])

# Heatmap
ha <- HeatmapAnnotation(Condition = meta.samples[sample,]$irAE,
                       Sample = sample,
                       col = list(Condition = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")),
                                  Sample = col.samples))

sample <- sample[order(match(sample, order.patient))]
pdf(paste0(dir.output, "Plots/DGE/TopClono/HeatmapAverage_Post_Top10_CD8_TEM_PB.pdf"), height = 9, width = 6)
Heatmap(data.matrix(mat2), show_column_names = F, top_annotation = ha,
        clustering_method_rows = "ward.D2",
        row_names_gp = gpar(fontsize = 6),
        column_order = names(sample))
dev.off()

# Select genes for OR (5c)

mat.sel <- small.seurat@assays$RNA@data
mat.sel <- mat.sel[setdiff(rownames(res), genes.to.remove),]
mat.sel <- sapply(unique(small.seurat$clonotype_sample), 
                  function(x) rowMeans(mat.sel[, rownames(small.seurat@meta.data[small.seurat$clonotype_sample == x,])]))

tab.samples <- table(sample)
samples.solo <- names(tab.samples[tab.samples == 1])
clono.solo <- names(sample[sample %in% samples.solo])
samples.multi <- setdiff(unique(sample), samples.solo)

mat.sel2 <- cbind(sapply(samples.multi, function(x) # patient agregation
  rowSums(mat.sel[, names(sample[sample == x])])), mat.sel[, clono.solo])
sel.genes <- intersect(rownames(mat.sel2[apply(mat.sel2, 1, function(x) sum(x > 0)) >= 3,]), #in at least 3 patients
                       rownames(mat.sel[apply(mat.sel, 1, function(x) sum(x >=0.2)) >= 8,])) #minimum value in at least 8 clonotypes (i.e. 5%)
length(sel.genes)