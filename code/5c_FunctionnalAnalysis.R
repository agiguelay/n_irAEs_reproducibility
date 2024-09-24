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

library(clusterProfiler)
library("org.Hs.eg.db", character.only = TRUE)
library(DOSE)
library(RColorBrewer)
library(ggplot2)

# Parameters --------------------------------------------------------------------
dir.output <- "..." 
col.samples <- ...
# Perform OR analysis for CD4 CTL ----------------------------------------------
res <- read.table(paste0(dir.output, "Tables/DGE_Top10_CD4_CTL_ControlvsirAE_Post_PB.txt"), header = T)
res <- rownames(res[res$FDR <= 0.01 & res$logFC <= -log2(1.5),])
sel <- read.table(paste0(dir.output, "Tables/GenesSel_Post_Top10_CD4_CTL_PB.txt"))$x
genes <- intersect(res, sel)

ora <- enrichGO(gene = genes,
               keyType = "SYMBOL",
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 1,
               qvalueCutoff  = 1,
               readable      = TRUE)
write.table(ora@result, paste0(dir.output, "Tables/ORA_GO_Top10_CD4_CTL_upirAE_PB.txt"), sep = '\t')

pdf(paste0(dir.output, "Plots/DGE/TopClono/DotPlot_GO_Top10_CD4_CTL_UpirAE_PB.pdf"),
    width = 9, height = 10)
print(dotplot(ora,
              showCategory = 20,
              orderBy = "qvalue",
              color = "qvalue") +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_line(color = "black", size = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", size = 0.8)) +
        scale_color_gradientn(colours = brewer.pal(9,"Spectral"))
)
dev.off()

# Perform OR analysis for CD8 TEM ----------------------------------------------
res <- read.table(paste0(dir.output, "Tables/DGE_Top10_CD8_TEM_ControlvsirAE_Post_PB.txt"), header = T)
res <- rownames(res[res$FDR <= 0.01 & res$logFC <= -log2(1.5),])
sel <- read.table(paste0(dir.output, "Tables/GenesSel_Post_Top10_CD8_TEM_PB.txt"))$x
genes <- intersect(res, sel)
ora <- enrichGO(gene = genes,
               keyType = "SYMBOL",
               OrgDb = org.Hs.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 1,
               qvalueCutoff  = 1,
               readable      = TRUE)
write.table(ora@result, paste0(dir.output, "Tables/ORA_GO_Top10_CD8_TEM_upirAE_PB.txt"), sep = '\t')

pdf(paste0(dir.output, "Plots/DGE/TopClono/DotPlot_GO_Top10_CD8_TEM_UpirAE_PB.pdf"),
    width = 9, height = 10)
print(dotplot(ora,
              showCategory = 20,
              orderBy = "qvalue",
              color = "qvalue") +
        theme(text = element_text(colour = "black", size = 15),
              axis.line = element_line(color = "black", size = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", size = 0.8)) +
        scale_color_gradientn(colours = brewer.pal(9,"Spectral"))
)
dev.off()