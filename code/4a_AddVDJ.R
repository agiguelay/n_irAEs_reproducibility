##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#              STEP4a: Add TCR information               #
#                                                        #
#    Gathers TCR data from all the samples               #
#    Defines the top clono for different compartments    #
#    Calculates the pool-normalized clono contribution   #
#    Calculates the proportion of top clono              #          
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
library(data.table)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(DescTools)
library(ggpubr)
library(robustrank)

# Functions --------------------------------------------------------------------

# @param df.tcr:
# @param seurat: seurat object
# @param name.cluster: name of the chosen compartment (for outputs)
# @param celltypes: vector with all the cell types included in the chosen compartment 
# @param name.sample: name of the metadata with the sample annotation
# @param name.annotation: name of the metadata with the cell type annotation
# @param dir.output: directory where all the outputs will be saved
# @param n1: defines the top n1 clones on which we focus (e.g. n1 = 10)
# @param n2: defines the top n2 clones to which we restrict the analysis (e.g. n2 = 50)
# @param min.exp: minimal number of cells where a TCR is detected to consider 
# the clone as expanded (e.g. min.exp = 2)
# Returns the seurat object with updated metadata + pool-normalized contribution table

TopClono <- function(df.tcr, seurat, name.cluster, celltypes, name.sample, 
                    name.annotation, dir.output, n1, n2, min.exp){
  #Initializes the future metadata (presence in the top n1; expanded or not)
  topclono <- expclono <- setNames(rep("No", ncol(seurat)), colnames(seurat))
  
  #Creates data.frames where each row corresponds to one of the top n2 clonotypes
  topn2.df <- data.frame(Rank = 1:n2)

  samples <- unique(seurat@meta.data[name.sample])[,1]
  
  for(sample in samples){
    #Subsets the cells (keeping only the ones present in @param celltypes and the respective sample)
    good.cells <- rownames(seurat@meta.data[seurat@meta.data[, name.sample] == sample & 
                                              seurat@meta.data[, name.annotation] %in% celltypes, ])
    small.seurat <- subset(seurat, cells = good.cells)
    inter <- intersect(unique(df.tcr$barcode), small.seurat$barcode)
    small.tcr <- df.tcr[df.tcr$barcode %in% inter,] #subset the tcr dataset
    
    # Counts the #times each clonotype is present and order on the frequency
    topn1_clonotype <- as.data.frame(table(small.tcr$clonotype_id))
    topn1_clonotype <- topn1_clonotype %>% arrange(desc(Freq))
    
    # Restricts to the top n2 clonotypes
    topn2 <- topn1_clonotype[1:n2,] 
    # Normalizes by the number of different clonotypes present (pool-normalized contribution)
    topn2 <- topn2$Freq/length(unique(topn1_clonotype$Var1))
    topn2.df <- cbind(topn2.df, topn2)
    
    top_n1_clonotype <- topn1_clonotype[1:n1,]
    exp.clono <- topn1_clonotype[topn1_clonotype$Freq >= min.exp,]
    
    # Updates the metadata and add them to the seurat object
    top.clono.cells <- rownames(small.tcr[small.tcr$clonotype_id %in% top_n1_clonotype$Var1,]) #cells corresponding to the top n1 clonotypes
    expanded.cells <- rownames(small.tcr[small.tcr$clonotype_id %in% exp.clono$Var1,]) #cells corresponding to expanded clonotypes
    
    expclono[expanded.cells] <- sample
    topclono[intersect(top.clono.cells, expanded.cells)] <- sample # in case the top clonotypes are actually not expanded
    
    seurat <- AddMetaData(seurat, topclono, col.name = paste0("Top", as.character(n1), "_", name.cluster))
    seurat <- AddMetaData(seurat, expclono, col.name = paste0("Exp_", name.cluster))
  }
  
  colnames(topn2.df)[-1] <- samples
  topn2.df[is.na(topn2.df)] <- 0
  
  fwrite(topn2.df, paste0(dir.output, "Tables/Top", as.character(n1),  "Clono_", name.cluster, ".txt"), sep = "\t")
  
  return(list(seurat, topn2.df))
}

# @param seurat: seurat object (might require subsetting)
# @param dir.output: directory where all the outputs will be saved
# @param analysis.ID: ID for the output files names
# @param name.cluster: name of the chosen compartment (for outputs)
# @param name.sample: name of the metadata with the sample annotation
# @param topn2.df: TopClono second output
# @param n1: defines the top n1 clones on which we focus (e.g. n1 = 10)
# @param n2: defines the top n2 clones to which we restrict the analysis (e.g. n2 = 50)
# @param color.by: name of the metadata to use to color the lines in the clonotype
# contribution plot
# @param linetype.by: name of the metadata to use to change the line type in the clonotype
# contribution plot
# @param boxplot.by: name of the metadata to use as a group for the AUC boxplot
# @param col.cond: vector of colors for the lines in the clonotype contribution plot
# @param col.samples: vector of colors for the different samples (AUC boxplot)
# Saves two plots (Clonotype's pool-normalized contribution plot and the respective AUC boxplot)

ClonoContribution <- function(seurat, dir.output, analysis.ID, name.cluster, 
                              name.sample, topn2.df, n1, n2, color.by, 
                              linetype.by, boxplot.by, col.cond, col.samples){
  
  # Prepares data 
  samples <- unique(seurat@meta.data[name.sample])[,1]
  keep.col <- colnames(topn2.df)[-1] %in% samples
  topn2.df <- topn2.df[, c(T, keep.col)]
  
  topn2.df.plot <- data.frame(Rank = rep(topn2.df$Rank, length(samples)),
                             Contribution = unlist(topn2.df[,-1]),
                             Sample = as.vector(sapply(samples, function(x) rep(x, n2)))
  )

  # Gives the color for the curves
  condition.col <- setNames(sapply(samples, function(x) seurat@meta.data[seurat@meta.data[name.sample] == x, color.by][1]), samples)
  condition.col <- sort(condition.col)
  
  # Gives the linetype for the curves
  condition.lt <- setNames(sapply(samples, function(x) seurat@meta.data[seurat@meta.data[name.sample] == x, linetype.by][1]), samples)
  condition.lt <- sort(condition.lt)
  
  # Gives the group for the boxplot
  condition.bp <- setNames(sapply(samples, function(x) seurat@meta.data[seurat@meta.data[name.sample] == x, boxplot.by][1]), samples)
  condition.bp <- sort(condition.bp)
  
  topn2.df.plot$Condition.col <- condition.col[topn2.df.plot$Sample]
  topn2.df.plot$Condition.line <- condition.lt[topn2.df.plot$Sample]
  
  # Curves 
  dir.create(paste0(dir.output, "Plots/VDJ/"), recursive = T)
  
  pdf(paste0(dir.output, "Plots/VDJ/ContriPoolnDistribution50_Top", as.character(n1),  "Clono_", name.cluster, "_", analysis.ID, ".pdf"), width = 8, height = 5)
  print(ggplot(topn2.df.plot, aes(x = as.numeric(Rank), y = Contribution, group = Sample)) +
          geom_line(linewidth = 0.8, aes(linetype = Condition.lt, color = Condition.col)) +
          xlab("Clonotype rank") +
          xlim(c(n2,1)) +
          scale_color_manual(values = col.cond) +
          scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
          theme(text = element_text(colour = "black", size = 15),
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  dev.off()

  # Calculate AUC 
  auc <- sapply(unique(topn2.df.plot$Sample), function(x) {
    df2b <- topn2.df.plot[topn2.df.plot$Sample == x,]
    return(AUC(df2b$Rank, df2b$Contribution))
  })
  
  auc.df <- data.frame(Sample = colnames(auc),
                      Pool = auc[2, ], #divided by the number of different TCR sequences
                      Condition.group = condition.bp[colnames(auc)],
                      Patient = sapply(colnames(auc), function(x) strsplit(x, "_")[[1]][1]))
  
  write.table(auc.df, paste0(dir.output, "Tables/Contri_Among50_forStats_AUC_", "Clono_", name.cluster, "_", analysis.ID, ".txt"))
  
  # Plot AUC
  pdf(paste0(dir.output, "Plots/VDJ/ContriPoolAUC50forpaper_Top", as.character(n1),  "Clono_", name.cluster, "_", analysis.ID, ".pdf"), width = 6, height = 5)
  print(ggplot(auc.df, aes(x = Condition.group, y = Pool)) +
          geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Sample), size = 3) +
          geom_boxplot(color = "black", linewidth = 1) +
          geom_point(size = 2) +
          theme_classic() +
          ylab("AUC") +
          scale_color_manual(values = col.samples) +
          stat_compare_means(method = "wilcox.test", paired = FALSE) +
          scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", linewidth = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", linewidth = 0.8)
          )
  )
  dev.off()

}


# @param seurat: seurat object
# @param dir.output: directory where all the outputs will be saved
# @param name.cluster: name of the chosen compartment (for outputs)
# @param name.sample: name of the metadata with the sample annotation
# @param name.annotation: name of the metadata with the cell type annotation
# @param name.cond: name of the metadata to group by (boxplot)
# @param n1: defines the top n1 clones on which we focus (e.g. n1 = 10)
# @param n2: defines the top n2 clones to which we restrict the analysis (e.g. n2 = 50)
# @param celltypes: vector with all the cell types included in the chosen compartment
# @param col.cond: vector of colors for the boxplot 
# Saves a serie of plots giving the proportion of cells in the top n1 clonotypes
# in the different cell types (boxplots)

CellTypeClono <- function(seurat, dir.output, analysis.ID, name.sample, name.cluster,
                          name.annotation, name.cond, n1, n2, celltypes, col.cond){

  # Subset the cells to the cell types of interest
  cells.clus <- rownames(seurat@meta.data[seurat@meta.data[, name.annotation] %in% celltypes, ])
  small.meta <- seurat@meta.data[cells.clus,]
  
  # Change the metadata (top n1 status) to something binary Y/N and 
  small.meta[, paste0("Top", as.character(n1), "_", name.cluster)] <- sapply(
    small.meta[, paste0("Top", as.character(n1), "_", name.cluster)], 
    function(x) ifelse(x == "No", yes = "No", no = "Yes"))

  # Counts the nb of cells in top n1 clones for each patient and each cell type
  ct.cl <- table(small.meta[[paste0("Top", as.character(n1), "_", name.cluster)]], 
                small.meta[[name.sample]], small.meta[[name.annotation]])
  
  # Serie of barplots (one for each cell type) representing the proportion of 
  # cells in the top n1
  pdf(paste0(dir.output, "Plots/VDJ/ExpandedClonesPerCellTypeProportion_Top", as.character(n1),  "_", name.cluster, "_", analysis.ID, ".pdf"), width = 10, height = 6)
  
  for(ct in celltypes){
    # Divides the counts by the total nb of cells in the patient for this specific cell type
    prop <- sweep(ct.cl[,, ct], 2, colSums(ct.cl[,, ct]), "/")
  
    # Prepares the data
    prop.plot <- data.frame(Sample = as.vector(sapply(colnames(prop), function(x) rep(x, nrow(prop)))),
                           Expanded = rep(rownames(prop), ncol(prop)),
                           Proportion = as.vector(prop))
    prop.plot$Condition <- condition[prop.plot$Sample]
    prop.plot$Condition <- factor(prop.plot$Condition, levels = names(col.cond))
    prop.plot$Patient <- sapply(prop.plot$Sample, function(x) strsplit(x, "_")[[1]][1])
    
    # Plot
    print(ggplot(prop.plot[prop.plot$Expanded == "Yes",], aes(x = Condition, y = Proportion, fill = Condition)) +
            geom_boxplot(color = "black", linewidth = 1) +
            geom_point(size = 2) +
            geom_line(aes(group = Patient)) +
            theme_classic() +
            ylab(paste0("Proportion of cells in top ", as.character(n1))) +
            scale_fill_manual(values = col.cond) +
            stat_compare_means(method = "wilcox.test", paired = FALSE) +
            scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
            theme(text = element_text(colour = "black", size = 15),
                  aspect.ratio = 1/1,
                  axis.line = element_line(color = "black", linewidth = 0.8),
                  axis.text = element_text(colour = "black", size = 15),
                  axis.ticks = element_line(color = "black", linewidth = 0.8))
    )
  }
  dev.off()
  
}

# Parameters -------------------------------------------------------------------
samples <- ... # vector with all the sample ID (have to be matching the cellranger output)
path.cr <- "..." # path for cellranger outputs
dir.output <- "..." 
n1 <- 10
n2 <- 50
min.exp <- 2
col.samples <- ...

# Prepare the TCR table from cellranger outputs --------------------------------

tcr.all = foreach(x = samples, .combine = "rbind") %do% {
  tcr.folder <- paste0(path.cr, x,"/outs/per_sample_outs/", x, "/vdj_t/")
  tcr <- fread(paste0(tcr.folder, "filtered_contig_annotations.csv"),
               header = T, stringsAsFactors = F, check.names = F, data.table = F)
  tcr <- tcr[order(tcr$cdr3),] #to be sure the same "consensus" are put assigned 
  # to consensus 1 or 2 across the cells of the same clonotype
  
  # Constructs four tables: 2 for the chains A and Beta (1rst consensus)
  # + 2 for the clonotypes possessing a second A or B chain (2nd consensus)
  tcra <- tcr[tcr$chain == "TRA", ]
  tcrb <- tcr[tcr$chain == "TRB", ]
  
  tcra1 <- tcra[!duplicated(tcra$barcode),]
  tcrb1 <- tcrb[!duplicated(tcrb$barcode),]
  tcra2 <- tcra[duplicated(tcra$barcode),]
  tcrb2 <- tcrb[duplicated(tcrb$barcode),]
  
  tcra <- left_join(tcra1, tcra2, by = c("barcode" = "barcode"), keep = F)
  tcrb <- left_join(tcrb1, tcrb2, by = c("barcode" = "barcode"), keep = F)
  
  # Only keep the barcode and clonotype columns. 
  # We'll get additional clonotype info from the clonotype table.
  tcra <- tcra[, c("barcode", "raw_clonotype_id.x","v_gene.x","d_gene.x","j_gene.x","v_gene.y","d_gene.y","j_gene.y")]
  colnames(tcra) <- c("barcode", "clonotype_id", "v_gene_A1", "d_gene_A1", "j_gene_A1", "v_gene_A2", "d_gene_A2", "j_gene_A2")
  
  tcrb <- tcrb[, c("barcode", "raw_clonotype_id.x","v_gene.x","d_gene.x","j_gene.x","v_gene.y","d_gene.y","j_gene.y")]
  colnames(tcrb) <- c("barcode", "clonotype_id", "v_gene_B1", "d_gene_B1", "j_gene_B1", "v_gene_B2", "d_gene_B2", "j_gene_B2")
  
  tcr <- full_join(tcra, tcrb, by = c("barcode" = "barcode", "clonotype_id" = "clonotype_id"))
  
  # Clonotype-centric info.
  clono <- fread(paste0(tcr.folder,"clonotypes.csv"), header = T, stringsAsFactors = F, check.names = F, data.table = F)
  
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
  
  rownames(tcr) <- paste0(as.character(x), "_", as.character(tcr$barcode))
  tcr$clonotype_sample <- paste0(as.character(tcr$clonotype_id), "_", as.character(x))
  tcr
}
write.table(tcr.all, paste0(dir.output, "Tables/TCR_Post_TRAandB.txt"), sep = "\t")

# Add it as metadata to the seurat object --------------------------------------
seurat.all <- readRDS(dir.output, "Objects/SeuratFinal_All.rds")

seurat.all <- AddMetaData(object = seurat.all, metadata = tcr.all)
seurat.all <- AddMetaData(object = seurat.all, metadata = 
                            sapply(seurat.all$Sample, function(x) strsplit(x, "_")[[1]][1]),
                          col.name = "Patient")

seurat.all$ClonoStatus <- sapply(seurat.all$clonotype_id, function(x) ifelse(is.na(x), 'Not available', "Available"))

# Project the TCR presence on a UMAP
pdf(paste0(dir.output, "Plots/VDJ/UMAP_Seurat_Integration_TCRavailable.pdf"), width = 8, height = 8)
print(DimPlot(seurat.all, reduction = "FinalUMAP", pt.size = 0.2, shuffle = T,
              group.by = "ClonoStatus", label = FALSE, raster = F) + NoLegend()  +
        scale_color_manual(values = setNames(c("#CCCC00","lavender"), c("Available", "Not available"))) +
        theme(text = element_text(colour = "black", size = 15),
              aspect.ratio = 1/1,
              axis.line = element_line(color = "black", linewidth = 0.8),
              axis.text = element_text(colour = "black", size = 15),
              axis.ticks = element_line(color = "black", linewidth = 0.8)
        ))

dev.off()


# Defines the top n1 clono of different compartments and plots the contribution ----
## CD8
res.topclono <- TopClono(df.tcr = tcr.all,
                        seurat = seurat.all,
                        name.cluster = "CD8",
                        celltypes = c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"),
                        name.sample = "Sample",
                        name.annotation = "FinalAnnotation",
                        dir.output = dir.output,
                        n1 = n1,
                        n2 = n2,
                        min.exp = min.exp) 

ClonoContribution(seurat = subset(res.topclono[[1]], subset = Pre_Post == "Post-ICI"),
                  topn2.df = res.topclono[[2]],
                  name.sample = "Sample",
                  name.cluster = "CD8",
                  dir.output = dir.output,
                  n1 = n1,
                  n2 = n2,
                  color.by = "irAE",
                  linetype.by = "irAE",
                  boxplot.by = "irAE",
                  analysis.ID = "Post",
                  col.cond = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")),
                  col.samples = col.samples)

## CD4
res.topclono <- TopClono(df.tcr = tcr.all,
                         seurat = res.topclono[[1]],
                         name.cluster = "CD4",
                         celltypes = c("CD4 Naive", "CD4 Proliferating", "CD4 CTL", "CD4 TEM", "CD4 TCM", "Treg"),
                         name.sample = "Sample",
                         name.annotation = "FinalAnnotation",
                         dir.output = dir.output,
                         n1 = n1,
                         n2 = n2,
                         min.exp = min.exp) 

ClonoContribution(seurat = subset(res.topclono[[1]], subset = Pre_Post == "Post-ICI"),
                  topn2.df = res.topclono[[2]],
                  name.sample = "Sample",
                  name.cluster = "CD4",
                  dir.output = dir.output,
                  n1 = n1,
                  n2 = n2,
                  color.by = "irAE",
                  linetype.by = "irAE",
                  boxplot.by = "irAE",
                  analysis.ID = "Post",
                  col.cond = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")),
                  col.samples = col.samples)

## CD8 TEM
res.topclono <- TopClono(df.tcr = tcr.all,
                         seurat = res.topclono[[1]],
                         name.cluster = "CD8_TEM",
                         celltypes = c("CD8 TEM"),
                         name.sample = "Sample",
                         name.annotation = "FinalAnnotation",
                         dir.output = dir.output,
                         n1 = n1,
                         n2 = n2,
                         min.exp = min.exp) 

## CD4 CTL 
res.topclono <- TopClono(df.tcr = tcr.all,
                         seurat = res.topclono[[1]],
                         name.cluster = "CD4_CTL",
                         celltypes = c("CDD CTL"),
                         name.sample = "Sample",
                         name.annotation = "FinalAnnotation",
                         dir.output = dir.output,
                         n1 = n1,
                         n2 = n2,
                         min.exp = min.exp) 


# Saves the results
saveRDS(res.topclono[[1]], paste0("Objects/SeuratFinal_All.rds"))

# Calculates the proportion of top n1 cells in each cell type ------------------
## CD8
CellTypeClono(seurat = subset(res.topclono[[1]], subset = Pre_Post == "Post-ICI"),
              dir.output = dir.output,
              analysis.ID = "Post",
              name.sample = "Sample",
              name.cluster = "CD8",
              name.annotation = "FinalAnnotation",
              name.cond = "irAE",
              n1 = n1,
              n2 = n2,
              celltypes = c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"),
              col.cond = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")))
## CD4
CellTypeClono(seurat = subset(res.topclono[[1]], subset = Pre_Post == "Post-ICI"),
              dir.output = dir.output,
              analysis.ID = "Post",
              name.sample = "Sample",
              name.cluster = "CD4",
              name.annotation = "FinalAnnotation",
              name.cond = "irAE",
              n1 = n1,
              n2 = n2,
              celltypes = c("CD4 Naive", "CD4 Proliferating", "CD4 CTL", "CD4 TEM", "CD4 TCM", "Treg"),
              col.cond = setNames(c("#CCCC00", "#CC0000"), c("Control", "irAE")))