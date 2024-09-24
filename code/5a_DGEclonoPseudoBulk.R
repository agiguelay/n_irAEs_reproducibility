##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#             STEP5a: Performs pseudobulk gene           #
#       differential expression analysis with edgeR      #
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################

# To perform on R 4.2.3 (edgeR 3.42.4)
rm(list = ls())

# Libraries --------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(ggrepel)
library(edgeR)

# Load functions ---------------------------------------------------------------

# @param seurat: seurat object
# @param dir.output: directory where all the outputs will be saved
# @param analysis.ID: ID for the output files names
# @param max.log.qval: force max value for the volcanot plot when the qvalue is too low
# @param max.qval: max value of qvalue to highlight a gene as differentially expressed
# @param min.fc: min value of fold change to highlight a gene as differentially expressed
# @param name.clono: name of the metadata with the unique clonotype annotation (ID + sample)
# @param name.cluster: name of the chosen compartment
# @param name.cond: name of the metadata used to oppose two groups
# @param col.cond: vector of colors for the two groups (VolcanoPlot)

DGE_Cond_PB = function(seurat, dir.output, analysis.ID, max.log.qval, max.qval, 
                       min.fc, name.clono, name.cluster, name.cond, col.cond){
  
  # Only keep the cells that are expanded/top 
  cells.to.keep <-  rownames(seurat@meta.data[seurat[[name.cluster]] != "No",])
  seurat <-  subset(seurat, cells = cells.to.keep)
  
  # Transforms single cell to pseudobulk
  y <- Seurat2PB(seurat, sample = name.cond, cluster = name.clono)
  
  # Filter out low detected genes
  keep.genes <- filterByExpr(y, group=y$samples$cluster, min.count = 5, min.total.count = 10) 
  y <- y[keep.genes, , keep = FALSE]
  
  # edgeR analysis
  y <- normLibSizes(y)
  condition <- factor(y$samples$sample)
  
  cl <- as.matrix(factor(y$samples$sample))
  cl[cl== levels(condition)[1],] = 1
  cl[cl==levels(condition)[2],] = 2
  c <- factor(cl)
  design <- model.matrix(~0+c)
  cm <- makeContrasts(c1-c2, levels = design)
  
  design <- model.matrix(~ 0 + condition)
  
  y <- estimateDisp(y, design, robust=TRUE)

  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, contrast=cm)
  out <- topTags(res, n = "Inf")$table

  write.table(out, paste0(dir.output, "Tables/DGE_", name.cluster, "_", levels(condition)[1], "vs", 
                          levels(condition)[2], "_", analysis.ID, "_PB.txt"), sep = "\t")

  df.plot <- data.frame(Gene = rownames(out),
                       log2FC = out$logFC,
                       logQval = -log10(out$FDR))
  
  if(sum(is.infinite(df.plot$logQval)) > 0){
    df.plot[is.infinite(df.plot$logQval),]$logQval <- max.log.qval
  }  
  
  df.plot$Status <- apply(df.plot, 1, function(x) 
    ifelse(test = as.numeric(x[3]) >= -log10(max.qval) & abs(as.numeric(x[2])) >= log2(min.fc),
           yes = ifelse(test = as.numeric(x[2]) >= log2(min.fc),
                        yes = paste0("Up in ", sort(unique(seurat@meta.data[[name.cond]]))[1]),
                        no = paste0("Up in ", sort(unique(seurat@meta.data[[name.cond]]))[2])),
           no = "No change"))
  
  # Removes TCR genes
  df.plot[grep("^TR", df.plot$Gene),]$Gene <- ""
  df.plot[df.plot$Status == "No change",]$Gene <- ""
  
  dir.create(paste0(dir.output, "Plots/DGE/"), recursive = T)
  
  pdf(paste0(dir.output, "Plots/DGE/VolcanoPlot_", levels(condition)[1], "vs", 
             levels(condition)[2], "_", name.cluster, "_", analysis.ID, "_PB.pdf"),
      width = 6, height = 4)
  print(ggplot(data = df.plot, mapping = aes(y = logQval, x = log2FC, color = Status, label = Gene)) +
          scale_color_manual(values = setNames(c("gray65", col.cond[1], col.cond[2]), 
                                               c("No change", paste0("Up in ", levels(condition)[1]), 
                                                 paste0("Up in ", levels(condition)[2]))))+
          geom_point(size = 0.8, alpha = 0.75, shape  = 19) +
          geom_text_repel(size = 2.5, max.overlaps = 50)+
          ylab("- log10(FDR)") +
          xlab("log2 Fold Change") +
          xlim(c(-max(abs(df.plot$log2FC))), max(abs(df.plot$log2FC)))+
          ylim(c(0, max(df.plot$logQval)+10))+
          theme_classic() +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8))
  )
  dev.off()
  
}

# Parameters --------------------------------------------------------------------
dir.output <- "..." 
samples <- ...
col.samples <- ...

# Load data --------------------------------------------------------------------
seurat.all <- readRDS(paste0(dir.output, "Objects/SeuratFinal_All.rds"))
seurat.all <- subset(seurat.all, subset = Pre_Post == "Post-ICI")

# DGE Top clones Ctrl vs irAE per clono ----------------------------------------

DGE_Cond_PB(seurat = seurat.all,
         dir.output = dir.output,
         analysis.ID = "Post",
         max.log.qval = 300, 
         max.qval = 0.01, 
         min.fc = 1.5, 
         name.cond = "irAE", 
         name.cluster = "Top10_CD8",
         name.clono = "clonotype_sample",
         col.cond = c("#CCCC00", "#CC0000"))

DGE_Cond_PB(seurat = seurat.all,
         dir.output = dir.output,
         analysis.ID = "Post",
         max.log.qval = 300, 
         max.qval = 0.01, 
         min.fc = 1.5, 
         name.cond = "irAE", 
         name.cluster = "Top10_CD4",
         name.clono = "clonotype_sample",
         col.cond = c("#CCCC00", "#CC0000"))

DGE_Cond_PB(seurat = seurat.all,
         dir.output = dir.output,
         analysis.ID = "Post",
         max.log.qval = 300, 
         max.qval = 0.01, 
         min.fc = 1.5, 
         name.cond = "irAE", 
         name.cluster = "Top10_CD4_CTL",
         name.clono = "clonotype_sample",
         col.cond = c("#CCCC00", "#CC0000"))

DGE_Cond_PB(seurat = seurat.all,
         dir.output = dir.output,
         analysis.ID = "Post",
         max.log.qval = 300, 
         max.qval = 0.01, 
         min.fc = 1.5, 
         name.cond = "irAE", 
         name.cluster = "Top10_CD8_TEM",
         name.clono = "clonotype_sample",
         col.cond = c("#CCCC00", "#CC0000"))
