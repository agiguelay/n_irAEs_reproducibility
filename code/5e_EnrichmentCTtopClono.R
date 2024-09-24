##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#          STEP5c: Over-representation analysis          #
#           of differentially expressed genes            #
#           (top clones) with clusterProfiler            #
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
library(dplyr)

# Parameters -------------------------------------------------------------------
dir.output <- "..."
pop.interest <- "CD4 CTL"
pop <- "Top10_CD4"

# Load data --------------------------------------------------------------------
seurat.all <- readRDS("/data/cephfs-1/work/projects/ludwig-ici-iraes/Results/Objects/SeuratFinal_All.rds")
seurat.all <- subset(seurat.all, subset = FinalAnnotation %in% c("CD4 CTL", "CD4 Proliferating", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg"))
# seurat.all <- subset(seurat.all, subset = FinalAnnotation %in% c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"))

nb.cells <- table(seurat.all$Sample, seurat.all$FinalAnnotation)
nb.cellsb <- cbind(nb.cells[, pop.interest], 
                  rowSums(nb.cells[, setdiff(colnames(nb.cells), pop.interest)]))

# Calculate enrichment pvalues -------------------------------------------------
seurat.all <- seurat.all[, seurat.all@meta.data[[pop]] != "No"]

res.enrichment = foreach(sample = unique(seurat.all$Sample), .combine = "rbind")%do%{
  small.seurat <- subset(seurat.all, subset = Sample == sample)
  nb.cells.clono <- table(small.seurat$clonotype_id, small.seurat$FinalAnnotation)
  
  if(!pop.interest %in% colnames(nb.cells.clono)){ #in case none of the cells are in pop.interest 
    nb.cells.clono <- cbind(nb.cells.clono, rep(0, nrow(nb.cells.clono)))
    colnames(nb.cells.clono)[ncol(nb.cells.clono)] = pop.interest
  }
  
  if(length(setdiff(colnames(nb.cells.clono), pop.interest)) > 1){ #normal case
    nb.cells.clonob <- cbind(nb.cells.clono[, pop.interest], rowSums(nb.cells.clono[, setdiff(colnames(nb.cells.clono), pop.interest)]))
  } else if (length(setdiff(colnames(nb.cells.clono), pop.interest)) == 1){ # in case all cells are in pop.interest or just one other celltype
    nb.cells.clonob <- cbind(nb.cells.clono[, pop.interest], nb.cells.clono[, setdiff(colnames(nb.cells.clono), pop.interest)])
  } else if (length(setdiff(colnames(nb.cells.clono), pop.interest)) == 0){ #in case all cells are in pop.interest
    nb.cells.clonob <- cbind(nb.cells.clono, rep(0, nrow(nb.cells.clono)))
  }
  topclono <- unique(small.seurat$clonotype_id)
  res <- sapply(topclono, function(clono) 1 - phyper(nb.cells.clonob[clono,1]-1, nb.cellsb[sample,1], nb.cellsb[sample,2], sum(nb.cells.clonob[clono,])))

  df <- data.frame(Sample = rep(sample, length(topclono)), 
                  Pvalue = res,
                  ClonoID = topclono)
}

res.enrichment$LogPvalue <- -log10(res.enrichment$Pvalue)

# Plot -------------------------------------------------------------------------
meta.samples <- read.table(paste0(dir.output, "../SamplesMetadata_All.txt"), sep = "\t", header = T)
meta.samples <- meta.samples[meta.samples$Dataset == "Maschmeyer",]

df.plot <- left_join(res.enrichment, meta.samples[, c(1, 6, 12)], by = "Sample")
df.plot$LogPvalue[is.infinite(df.plot$LogPvalue)] <- 16
df.plot <- df.plot[df.plot$Pre_Post == "Post-ICI",]
df.plot %>% group_by(irAE) %>% 
  summarize(Mean = mean(LogPvalue), Median = median(LogPvalue))

pdf(paste0(dir.output, "Plots/VDJ/BoxPlot_NbClonoEnriched_Exp_", pop.interest ,".pdf"), width = 6, height = 6)
print(ggplot(df.plot, aes(x = irAE, y = LogPvalue, fill = irAE)) +
        geom_boxplot(color = "black", linewidth = 1, outlier.shape = NA)  +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.2), dotsize = 1, stackratio = 0.1) +
        scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        theme_classic() +
        ylab("-log10(p.value)") +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)

print(ggplot(df.plot, aes(x = irAE, y = LogPvalue, group.by = irAE)) +
        geom_boxplot(color = "black", linewidth = 1, outlier.shape = NA, fill = "white")  +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.2), dotsize = 1, stackratio = 0.1, aes(fill = Sample)) +
        scale_fill_manual(values = col.samples) +
        theme_classic() +
        ylab("-log10(p.value)") +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8))
)
dev.off()

