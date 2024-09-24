##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#       STEP3a: Compare cell type compositions           #
#                                                        #
#    Compares cell type cell type compositions between   #
#    ctrl/n-irae groups or pre/post treatment for the    # 
#    ctrl group                                          #
#    Calculates a B cell maturity score                  #
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
library(ggpubr)
library(dplyr)
library(robustrank)

# Load data --------------------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))

meta.samples <- read.table(dir.output, "SamplesMetadata_All.txt", sep = "\t", header = T)
meta.samples <- meta.samples[meta.samples$Dataset == "Maschmeyer",]

col.palette <- c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")

meta.samples <- meta.samples[order(meta.samples$irAE),] #To have consistent colors between controls and irAE patients
col.samples <- setNames(col.palette[1:nrow(meta.samples)], meta.samples$Sample)
col.celltypes <- setNames(col.palette[1:length(unique(seurat.all$predicted.celltype.l2))], sort(unique(seurat.all$predicted.celltype.l2)))


# Cell type proportion ---------------------------------------------------------
## Between Ctrl Post and irAE Post----
nb.cells <- table(seurat.all$FinalAnnotation, seurat.all$orig.ident, seurat.all$Pre_Post)
nb.post <- nb.cells[,, "Post-ICI"]
nb.post <- nb.post[, colSums(nb.post) > 0] # Select only Post samples
prop.post <- sweep(nb.post, 2, STATS = colSums(nb.post), FUN = "/")
prop.post <- data.frame(Proportion = as.vector(prop.post), CellType = rep(rownames(prop.post), ncol(prop.post)), 
                       Patient = as.vector(sapply(colnames(prop.post), function(x) rep(x, nrow(prop.post)))))
prop.post2 <- left_join(prop.post, meta.samples, by = c("Patient" = "Sample"), keep = F)

dir.create(paste0(dir.output, "Plots/Proportions/"), recursive = T)

pdf(paste0(dir.output, "Plots/Proportions/CellTypeProportionSingle.pdf"), width = 6, height = 6)
for(ct in unique(prop.post2$CellType)){
  small.prop.post2 = prop.post2[which(prop.post2$CellType == ct),]
  print(ggplot(small.prop.post2, aes(x = irAE, y = Proportion, fill = irAE)) +
          geom_boxplot(color = "black", linewidth = 1) + ggtitle(ct) +
          geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  
}
dev.off()

# B cells ----------------------------------------------------------------------
b.cells <- table(seurat.all$FinalAnnotation, seurat.all$orig.ident)
b.cells <- b.cells[c("B intermediate", "B memory", "B naive", "Plasmablast"),]
b.prop <- sweep(b.cells, 2, colSums(b.cells), "/")

df.plot.1 <- data.frame(Proportion = as.vector(b.prop), 
                       CellType = rep(rownames(b.prop), ncol(b.prop)),
                       Sample = as.vector(sapply(colnames(b.prop), function(x) rep(x, nrow(b.prop)))))
df.plot.1 <- left_join(df.plot.1, meta.samples, by = c("Sample" = "Sample"), keep = F)

df.plot.1$CellType <- factor(df.plot.1$CellType, levels = c("Plasmablast", "B memory", "B intermediate", "B naive"))
df.plot.1$Sample <- factor(df.plot.1$Sample, levels = unique(meta.samples$Sample))
df.plot.1$Patient <- sapply(as.character(df.plot.1$Sample), function(x) strsplit(x, "_")[[1]][1])

pdf(paste0(dir.output, "Plots/Proportions/CellTypeProportionBcells_Single.pdf"), width = 6, height = 6)
for(ct in unique(df.plot.1$CellType)){
  small.df.plot.1 <- df.plot.1[which(df.plot.1$CellType == ct & df.plot.1$Pre_Post == "Post-ICI"),]
  print(ggplot(small.df.plot.1, aes(x = irAE, y = Proportion, fill = irAE)) +
          geom_boxplot(color = "black", linewidth = 1) + ggtitle(ct) +
          geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  
}
dev.off()

# Calculates a b cell maturity score
b.score <- apply(b.prop, 2, function(x) 0.5*x[1] + 1*x[2] + 1*x[4])
df.plot.2 <- data.frame(Sample = names(b.score), Score = b.score)
df.plot.2 <- left_join(df.plot.2, meta.samples, by = c("Sample" = "Sample"), keep = F)
df.plot.2$Sample <- factor(df.plot.2$Sample, levels = unique(meta.samples$Sample))
rownames(df.plot.2) <- df.plot.2$Sample
df.plot.2$Patient <- sapply(as.character(df.plot.2$Sample), function(x) strsplit(x, "_")[[1]][1])

pvalue2 <- pm.wilcox.test(df.plot.2[x.pair, "Score"], df.plot.2[y.pair, "Score"], df.plot.2[x.plus, "Score"], df.plot.2[y.plus, "Score"])$p.value

pdf(paste0(dir.output, "Plots/Proportions/BcellScore.pdf"), width = 5, height = 4)
print(ggplot(df.plot.2[df.plot.2$Pre_Post == "Post-ICI",], aes(x = irAE, y = Score, fill = irAE)) +
        geom_boxplot(color = "black", linewidth = 1) + ggtitle("B cell maturity score") +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
        scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
        theme_classic() +
        stat_compare_means(method = "wilcox.test", paired = FALSE) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        )
)

print(ggplot(df.plot.2[df.plot.2$irAE == "Control",], aes(x = Pre_Post, y = Score, fill = Pre_Post)) +
        geom_boxplot(color = "black", linewidth = 1) + ggtitle("B cell maturity score") +
        geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
        scale_fill_manual(values = setNames(c("#333399", "#CC9900"), c("Pre-ICI", "Post-ICI"))) +
        theme_classic() +
        geom_line(aes(group = Patient)) +
        xlim(c("Pre-ICI", "Post-ICI")) +
        annotate("text", x = 1.3, y = 0.8, label = paste0("part. paired Wilcoxon, p = ", as.character(round(pvalue2, digits = 2)))) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        )
)
dev.off()

# CD4 --------------------------------------------------------------------------
cd4.cells <- table(seurat.all$FinalAnnotation, seurat.all$orig.ident)
cd4.cells <- cd4.cells[c("CD4 Naive", "CD4 CTL", "CD4 TCM", "CD4 TEM", "Treg", "CD4 Proliferating"),]
cd4.prop <- sweep(cd4.cells, 2, colSums(cd4.cells), "/")

df.plot.3 <- data.frame(Proportion = as.vector(cd4.prop), 
                       CellType = rep(rownames(cd4.prop), ncol(cd4.prop)),
                       Sample = as.vector(sapply(colnames(cd4.prop), function(x) rep(x, nrow(cd4.prop)))))
df.plot.3 <- left_join(df.plot.3, meta.samples, by = c("Sample" = "Sample"), keep = F)

df.plot.3$CellType <- factor(df.plot.3$CellType, levels = c("CD4 Naive", "CD4 CTL", "CD4 TCM", "CD4 TEM", "Treg", "CD4 Proliferating"))
df.plot.3$Sample <- factor(df.plot.3$Sample, levels = unique(meta.samples$Sample))
df.plot.3$Patient <- sapply(as.character(df.plot.3$Sample), function(x) strsplit(x, "_")[[1]][1])

pdf(paste0(dir.output, "Plots/Proportions/CD4_Proportion.pdf"), width = 8, height = 8)
print(ggplot(df.plot.3[df.plot.3$Pre_Post == "Post-ICI",], aes(x = Patient, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", color = "black", linewidth = 1) +
        theme_classic() +
        scale_fill_manual(values = col.celltypes) +
        coord_flip() +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_blank(),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_blank()
        ) 
)

print(ggplot(df.plot.3[df.plot.3$irAE == "Control",], aes(x = Patient, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", color = "black", linewidth = 1) +
        theme_classic() +
        scale_fill_manual(values = col.celltypes) +
        coord_flip() +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_blank(),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_blank()
        ) 
)
dev.off()

pdf(paste0(dir.output, "Plots/Proportions/CellTypeProportionCD4_Single.pdf"), width = 6, height = 6)
for(ct in unique(df.plot.3$CellType)){
  small.df.plot.3 <- df.plot.3[which(df.plot.3$CellType == ct & df.plot.3$Pre_Post == "Post-ICI"),]
  print(ggplot(small.df.plot.3, aes(x = irAE, y = Proportion, fill = irAE)) +
          geom_boxplot(color = "black", linewidth = 1) + ggtitle(ct) +
          geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  
}
dev.off()

# CD8 --------------------------------------------------------------------------
cd8.cells <- table(seurat.all$FinalAnnotation, seurat.all$orig.ident)
cd8.cells <- cd8.cells[c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating"),]
cd8.prop <- sweep(cd8.cells, 2, colSums(cd8.cells), "/")

df.plot.4 <- data.frame(Proportion = as.vector(cd8.prop), 
                       CellType = rep(rownames(cd8.prop), ncol(cd8.prop)),
                       Sample = as.vector(sapply(colnames(cd8.prop), function(x) rep(x, nrow(cd8.prop)))))
df.plot.4 <- left_join(df.plot.4, meta.samples, by = c("Sample" = "Sample"), keep = F)

df.plot.4$CellType <- factor(df.plot.4$CellType, levels = c("CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating"))
df.plot.4$Sample <- factor(df.plot.4$Sample, levels = unique(meta.samples$Sample))
df.plot.4$Patient <- sapply(as.character(df.plot.4$Sample), function(x) strsplit(x, "_")[[1]][1])

pdf(paste0(dir.output, "Plots/Proportions/CD8_Proportion.pdf"), width = 8, height = 8)
print(ggplot(df.plot.4[df.plot.4$Pre_Post == "Post-ICI",], aes(x = Patient, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", color = "black", linewidth = 1) +
        theme_classic() +
        scale_fill_manual(values = col.celltypes) +
        coord_flip() +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_blank(),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_blank()
        ) 
)

print(ggplot(df.plot.4[df.plot.4$irAE == "Control",], aes(x = Patient, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", color = "black", linewidth = 1) +
        theme_classic() +
        scale_fill_manual(values = col.celltypes) +
        coord_flip() +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_blank(),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_blank()
        ) 
)
dev.off()

pdf(paste0(dir.output, "Plots/Proportions/CellTypeProportionCD8_Single.pdf"), width = 6, height = 6)
for(ct in unique(df.plot.4$CellType)){
  small.df.plot.4 <- df.plot.4[which(df.plot.4$CellType == ct & df.plot.4$Pre_Post == "Post-ICI"),]
  print(ggplot(small.df.plot.4, aes(x = irAE, y = Proportion, fill = irAE)) +
          geom_boxplot(color = "black", linewidth = 1) + ggtitle(ct) +
          geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  
}
dev.off()

# Proportion populations -------------------------------------------------------
cells <- table(seurat.all$FinalAnnotation, seurat.all$orig.ident)
cells <- sweep(cells, 2, colSums(cells), "/")
df.plot.7 <- data.frame(Proportion = c(colSums(cells[c("B intermediate", "B memory", "B naive", "Plasmablast"),]), colSums(cells[c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"),]), 
                                      colSums(cells[c("CD4 CTL", "CD4 Proliferating", "CD4 Naive", "CD4 TCM", "CD4 TEM", "Treg"),]), colSums(cells[c("NK", "NK Proliferating", "NK_CD56bright"),]), 
                                      colSums(cells[c("CD16 Mono", "CD14 Mono"),])),
                       Sample = rep(colnames(b.prop), 5),
                       CellType = as.vector(sapply(c("B cells", "CD8 T cells", "CD4 T cells", "NK cells", "Monocytes"), function(x) rep(x, ncol(b.prop)))))
df.plot.7$Patient <- sapply(df.plot.7$Sample, function(x) strsplit(x, "_")[[1]][1])
df.plot.7 <- left_join(df.plot.7, meta.samples, by = c("Sample" = "Sample"), keep = F)

pdf(paste0(dir.output, "Plots/Proportions/CellTypeProportionBigPolulations_Single.pdf"), width = 6, height = 6)
for(ct in unique(df.plot.7$CellType)){
  small.df.plot.7 <- df.plot.7[which(df.plot.7$CellType == ct & df.plot.7$Pre_Post == "Post-ICI"),]
  print(ggplot(small.df.plot.7, aes(x = irAE, y = Proportion, fill = irAE)) +
          geom_boxplot(color = "black", linewidth = 1) + ggtitle(ct) +
          geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.75), dotsize = 1) +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15, angle = 45, hjust = 1),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  
}
dev.off()