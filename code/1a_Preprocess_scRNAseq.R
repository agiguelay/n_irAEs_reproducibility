##########################################################
#                                                        #
#       Neurological immune-related adverse effects      #
#         induced by Immuno Checkpoint Inhibition        #
#                                                        #
##########################################################
#                                                        #
#        STEP1a: Quality control of scRNAseq data        #
#                                                        #
#    Creates Seurat objects for each cell ranger output  #
#    Plots QC for individual samples                     #
#    Plots QC across all the samples                     #            
#                                                        #
##########################################################
#                                                        #
#                    Ambre Giguelay                      #
#                                                        #
##########################################################


rm(list = ls())

# Libraries --------------------------------------------------------------------

library(foreach)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(data.table)
library(ggplot2)

# Functions --------------------------------------------------------------------

# @param data: path of the cellranger sample_filtered_feature_bc_matrix output
# @param sample.ID: sample's name/ID
# @param dir.output: directory where all the outputs will be saved
# saves QC plot and returns a table with metadata

QC_scRNAseq_per_sample <- function(data, sample.ID, dir.output, ...){
  
  cat(paste0("Preprocessing of library ", sample.ID, "\n"))
  data.seurat <-  Read10X(data.dir = data)
  
  cat(paste0("Dataset before filtering: ", "\n", "\t", "- ", nrow(data.seurat), " genes", "\n", "\t", "- ", ncol(data.seurat), " cells", "\n"))
  data.seurat <-  CreateSeuratObject(counts = data.seurat, project = sample.ID)
  data.seurat[["percent.mt"]] <-  PercentageFeatureSet(data.seurat, pattern = "^MT-")
  data.seurat[["percent.rbs"]] <-  PercentageFeatureSet(data.seurat, pattern = "^RPS|^RPL")
  
  dir.create(paste0(dir.output, "Objects/"), recursive = T)
  saveRDS(object = data.seurat, file = paste0(dir.output, "Objects/SeuratObject_", sample.ID, ".rds"))
  
  df <-  data.frame(data.seurat@meta.data)
  
  dir.create(paste0(dir.output, "Plots/scRNAseq_QC/"), recursive = T)
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nFeatures_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = nFeature_RNA, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#CC9933", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# detected features (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nUMIs_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = nCount_RNA, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#CCCC33", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# UMIs (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_MTC_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = percent.mt, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#003366", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("% mitochondrial content") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_RBS_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = percent.rbs, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#339999", color = "black") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("% ribosomal content") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15),
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_nFeatures_vs_nUMIs_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = nFeature_RNA, x = nCount_RNA)) +
          geom_point(size = 0.8, color = "#CC9933", alpha = 0.5, shape  = 19) +
          ylab("# detected features (RNA)") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_MTC_vs_nUMIs_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = percent.mt, x = nCount_RNA)) +
          geom_point(size = 0.8, color = "#003366", alpha = 0.5, shape  = 19) +
          ylab("% mitochondrial content") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_RBS_vs_nUMIs_", sample.ID, ".pdf"), width = 4, height = 4)
  print(ggplot(data = df, mapping = aes(y = percent.rbs, x = nCount_RNA)) +
          geom_point(size = 0.8, color = "#339999", alpha = 0.5, shape  = 19) +
          ylab("% ribosomal content") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                aspect.ratio = 1/1,
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  return(df)
}



# @param data.all: dataframe with the metadata of all individual samples
# @param dir.output: directory where all the outputs will be saved
# @param max.mtc: maximal value for mitochondrial content in a cell
# @param max.rbs: maximal value for ribosomal content in a cell
# @param min.features: minimal number of features detected in a cell
# @param min.umis: minimal number of UMIs detected in a cell
# @param max.umis: maximal number of UMIs detected in a cell (use a high number to avoid messing up with doublet removal)
# saves QC plots

QC_scRNAseq_across_sample = function(data.all, dir.output, max.mtc, max.rbs, min.features, min.umis, max.umis){
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nFeatures_AllSamples.pdf"), width = 30, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = nFeature_RNA, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#CC9933", color = "black") +
          geom_hline(yintercept = min.features, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 3, fill = "snow", shape = 23) +
          ylab("# detected features (RNA)") +
          xlab("Library") +
          theme(text = element_text(colour = "black", size = 13), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 13, angle = 30, hjust = 1),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  df.all$Data = "scRNAseq" 
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nFeatures_AllSamplesCombined.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = nFeature_RNA, x = Data)) +
          geom_violin(size = 0.8, fill = "#CC9933", color = "black") +
          geom_hline(yintercept = min.features, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# detected features (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nUMIs_AllSamples.pdf"), width = 30, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = nCount_RNA, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#CCCC33", color = "black") +
          geom_hline(yintercept = min.umis, size = 0.8, color = "black", linetype = "dashed") +
          # geom_hline(yintercept = max.umis, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 3, fill = "snow", shape = 23) +
          ylab("# UMIs (RNA)") +
          xlab("Library") +
          theme(text = element_text(colour = "black", size = 13), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 13, angle = 30, hjust = 1),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_nUMIs_AllSamplesCombined.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = nCount_RNA, x = Data)) +
          geom_violin(size = 0.8, fill = "#CCCC33", color = "black") +
          geom_hline(yintercept = min.umis, size = 0.8, color = "black", linetype = "dashed") +
          # geom_hline(yintercept = max.umis, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("# UMIs (RNA)") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_MTC_AllSamples.pdf"), width = 30, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.mt, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#003366", color = "black") +
          geom_hline(yintercept = max.mtc, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 3, fill = "snow", shape = 23) +
          ylab("% mitochondrial content") +
          xlab("Library") +
          theme(text = element_text(colour = "black", size = 13), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 13, angle = 30, hjust = 1),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_MTC_AllSamplesCombined.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.mt, x = Data)) +
          geom_violin(size = 0.8, fill = "#003366", color = "black") +
          geom_hline(yintercept = max.mtc, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("% mitochondrial content") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_MRBS_AllSamples.pdf"), width = 30, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.rbs, x = orig.ident)) +
          geom_violin(size = 0.8, fill = "#339999", color = "black") +
          geom_hline(yintercept = max.rbs, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 3, fill = "snow", shape = 23) +
          ylab("% ribosomal content") +
          xlab("Library") +
          theme(text = element_text(colour = "black", size = 13), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 13, angle = 30, hjust = 1),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Vln_RBS_AllSamplesCombined.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.rbs, x = Data)) +
          geom_violin(size = 0.8, fill = "#339999", color = "black") +
          geom_hline(yintercept = max.rbs, size = 0.8, color = "black", linetype = "dashed") +
          stat_summary(fun.y = median, geom = "point", size = 5, fill = "snow", shape = 23) +
          ylab("% ribosomal content") +
          xlab("") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_nFeatures_vs_nUMIs_AllSamples.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = nFeature_RNA, x = nCount_RNA)) +
          geom_point(size = 0.5, color = "#CC9933", alpha = 0.5, shape  = 19) +
          geom_vline(xintercept = min.umis, size = 0.8, color = "black", linetype = "dashed") +
          # geom_vline(xintercept = max.umis, size = 0.8, color = "black", linetype = "dashed") +
          geom_hline(yintercept = min.features, size = 0.8, color = "black", linetype = "dashed") +
          ylab("# detected features (RNA)") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_MTC_vs_nUMIs_AllSamples.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.mt, x = nCount_RNA)) +
          geom_point(size = 0.5, color = "#003366", alpha = 0.5, shape  = 19) +
          geom_vline(xintercept = min.umis, size = 0.8, color = "black", linetype = "dashed") +
          # geom_vline(xintercept = max.umis, size = 0.8, color = "black", linetype = "dashed") +
          geom_hline(yintercept = max.mtc, size = 0.8, color = "black", linetype = "dashed") +
          ylab("% mitochondrial content") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()
  
  pdf(paste0(dir.output, "Plots/scRNAseq_QC/Scatter_RBS_vs_nUMIs_AllSamples.pdf"), width = 4, height = 4)
  print(ggplot(data = df.all, mapping = aes(y = percent.rbs, x = nCount_RNA)) +
          geom_point(size = 0.5, color = "#339999", alpha = 0.5, shape  = 19) +
          geom_vline(xintercept = min.umis, size = 0.8, color = "black", linetype = "dashed") +
          # geom_vline(xintercept = max.umis, size = 0.8, color = "black", linetype = "dashed") +
          geom_hline(yintercept = max.mtc, size = 0.8, color = "black", linetype = "dashed") +
          ylab("% ribosomal content") +
          xlab("# UMIs (RNA)") +
          theme(text = element_text(colour = "black", size = 15), 
                axis.line = element_line(color = "black", size = 0.8),
                axis.text = element_text(colour = "black", size = 15),
                axis.ticks = element_line(color = "black", size = 0.8)))
  dev.off()

}
