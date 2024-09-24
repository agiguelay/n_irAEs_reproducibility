##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#     STEP3e: Calculate AUCell scores for exhaustion,    #
#               naivety and cytotoxicity                 #
#                                                        #
#    For CD4 T cells or CD8 T cells:                     #
#    Calculates signature score with AUCells             #
#    Diverse plots                                       #
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
library(AUCell)
library(ggpubr)
library(foreach)

# Parameters -------------------------------------------------------------------
dir.output <- "..."
ct <- "CD8" #or CD4

col.palette <- c("#333399", "#6666CC", "#9999FF", "#CCCCFF", "#e3b3e3", "#b37bb3", "#993366", "#82032d", "#bd2f00", "#e05b02", "#FF9900", "#e8ba02", 
                "#b5b500", "#669900", "#52750c", "#009966", "cyan4", "turquoise3", "paleturquoise3", "powderblue", "lightskyblue3", "deepskyblue4", 
                "lightsteelblue4", "#FFCC99", "#CC9966", "#996633", "#663300", "#CC9999", "#996666", "#FF9999", "#CC6666", "black")

col.ct <- setNames(col.palette[12:15], c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"))
# col.ct <- setNames(col.palette[7:11], c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM"))
size.win <- 100

# Load data =-------------------------------------------------------------------
seurat <- readRDS(paste0(dir.output, "Objects/SeuratIntegrated_", ct, ".rds"))
pt <- seurat$Pseudotime
umap <- Embeddings(seurat, reduction = "umap")

seurat <- readRDS("/data/cephfs-1/work/projects/ludwig-ici-iraes/Results/Objects/SeuratFinal_All.rds") #has the complete metadata
seurat <- subset(seurat, cells = names(pt))
seurat <- AddMetaData(seurat, pt, col.name = "Pseudotime")

exprMatrix <- seurat@assays$RNA@counts

# Gene sets --------------------------------------------------------------------
tex <- c("CTLA4", "TIGIT", "PD1", "LAG3", "TIM3")
cyto <- c("NKG7", "CCL4", "CST7", "PRF1", "GZMA", "GZMB", "IFNG", "CCL3", "FGFBP2")# https://doi.org/10.1002/advs.202101447
naive <- c("CCR7", "SELL", "LEF1", "SELL")# https://doi.org/10.1002/advs.202101447

geneSets <- list("tex" = tex, "cyto" = cyto, "naive" = naive)

# Runs AUCell -------------------------------------------------------------------
cells_AUC <- AUCell_run(exprMatrix, geneSets)
res <- getAUC(cells_AUC)
write.table(res, paste0(dir.output, "Tables/AUCell_", ct, ".txt"), sep = "\t")

# Plots ------------------------------------------------------------------------

## umap -----

df.plot <- cbind(umap, as.data.frame(t(data.matrix(res))))
df.plot$Pseudotime <- pt[rownames(df.plot)]
df.plot$irAE <- seurat$irAE[rownames(df.plot)]
df.plot$CellType <- seurat$FinalAnnotation[rownames(df.plot)]

dir.create(paste0(dir.output, "Plots/GeneProgram/"), recursive = T)

 for(sig in names(geneSets)){
   df.plot$Signature <- df.plot[, sig]
   pdf(paste0(dir.output, "Plots/GeneProgram/UMAP_Score_", sig, "_", ct, ".pdf"), width = 5, height = 5)
   print(ggplot(df.plot, aes(x = UMAP_1, y = UMAP_2, color = Signature)) +
           geom_point(size = 0.3) +
           scale_color_gradientn(colors = c("#333399", "#0099CC", "#99CCFF", "#FFFFCC", "#FFCC33", "#FF9900", "#CC3300")) +
           theme_classic() +
           theme(text = element_text(colour = "black", size = 15),
                 aspect.ratio = 1/1,
                 axis.line = element_line(color = "black", linewidth = 0.8),
                 axis.text = element_text(colour = "black", size = 15),
                 axis.ticks = element_line(color = "black", linewidth = 0.8)
           ))
   
   dev.off()
 }


## scatter plots ----
pairs <- combn(names(geneSets), 2)

apply(pairs, 2, function(x) {
  sig1 <- x[1]
  sig2 <- x[2]
  df.plot$Signature1 <- df.plot[, sig1]
  df.plot$Signature2 <- df.plot[, sig2]
  pdf(paste0(dir.output, "Plots/GeneProgram/ScatterPlot_Score_", sig1, "vs", sig2,  "_", ct, ".pdf"), width = 5, height = 5)
  print(ggplot(df.plot, aes(x = Signature1, y = Signature2, color = irAE)) +
          geom_point(alpha = 0.5) +
          scale_color_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          xlab(sig1) + ylab(sig2)+
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          ))
  
  print(ggplot(df.plot, aes(x = Signature1, y = Signature2, color = CellType)) +
          geom_point(alpha = 0.5) +
          scale_color_manual(values = col.ct) +
          theme_classic() +
          xlab(sig1) + ylab(sig2)+
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          ))
  dev.off()

})

## pseudotime ----
df.plot2 <- df.plot[!is.na(df.plot$Pseudotime),]
df.plot2 <- df.plot2[order(df.plot2$Pseudotime),]

df.plot.all = foreach(sig = names(geneSets), .combine = "rbind") %do% {
  df.plot2$Signature <- df.plot2[, sig]
  df.plot3 <- df.plot2[size.win:(nrow(df.plot2)-size.win),]
  df.plot3$SignatureMean <- sapply(size.win:(nrow(df.plot2)-size.win), function(x) mean(df.plot2$Signature[(x-size.win):(x+size.win)]))
  df.plot3$Score <- rep(sig, nrow(df.plot3))
  df.plot3
}

pdf(paste0(dir.output, "Plots/GeneProgram/Pseudotime_Score_Allscores_", ct, ".pdf"), width = 15, height = 5)
print(ggplot(df.plot.all, aes(x = Pseudotime, y = SignatureMean, color = Score)) +
        geom_line(linewidth = 1) +
        scale_color_manual(values = setNames(c("#CCCCFF", "#FF9933", "#FFCC99"), c("naive1", "cyto1", "tex3"))) +
        theme_classic() +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/2,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))

dev.off()

## violin plots -----

for(score in c("tex", "cyto", "naive")){
  df.plot$Score <- df.plot[, score]
  
  pdf(paste0(dir.output, "Plots/GeneProgram/ViolinPlot_Score_", score, "_", ct, ".pdf"), width = 15, height = 5)
  print(ggplot(df.plot, aes(x = CellType, y = Score, fill = irAE)) +
          geom_violin() +
          scale_fill_manual(values = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE"))) +
          theme_classic() +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/2,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          ))
  
  dev.off()
  
}

## heatmap -----

pdf(paste0(dir.output, "Plots/GeneProgram/Heatmap_Scores_", ct, "_test.pdf"), width = 8, height = 2.2)
for(score in c("Pseudotime", "cyto", "tex", "naive")){
  print(score)
  df.plot$Score <- df.plot[, score]
  stats.df <- df.plot %>% group_by(CellType, irAE) %>% summarize(Mean = mean(Score, na.rm = T), Median = median(Score, na.rm = T))
  print(stats.df)

  median.df <- stats.df$Median
  mean.df <- stats.df$Mean
  
  df <- rbind(median.df, rep(min(mean.df), length(median.df)))
  colnames(df) <- paste0(stats.df$CellType, "_", stats.df$irAE)
  rownames(df) <- c(score, "fake")
  # df <- df[, c(3:8, 1:2)]
  
  print(pheatmap(df,
                 cluster_rows = F,
                 cluster_cols = F,
                 gaps_row = 1,
                 gaps_col = seq(2, length(median.df), 2),
                 show_rownames = F,
                 show_colnames = T, 
                 display_numbers = T,
                 number_format = "%.2f",
                 number_color = "white",
                 fontsize_number = 12,
                 border_color = "white",
                 cell_height = 0.5))

}
dev.off()



