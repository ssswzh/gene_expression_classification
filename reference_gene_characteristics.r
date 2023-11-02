#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen

library(ggplot2, quietly=TRUE)
theme_set(theme_classic())
library(ggpubr)


if ( FALSE ) {
  setwd("F://项目//TEP子宫肉瘤//tpm_matrix//merge_batch1_batch2//test_reference_genes")
  gene_file <- "reference_gene.list"
  matrix_file <- "..//test_limma//count.limma_normalized.tsv"
  sample_group_file <- "..//conditions.with_cofactor.all_sample_two_batch.txt"
}


ReadFiles <- function ( matrix_file, sample_group_file, gene_file=NULL) {
  norm_matrix <- read.delim(matrix_file, header=T, sep="\t", row.names=1, stringsAsFactors=F, check.names=F)
  sample_group <- read.delim(sample_group_file, header=T, row.names=1, sep="\t", stringsAsFactors=F, check.names=F)
  # only include samples in sample file and reorder
  shared_id <- sort(intersect(rownames(sample_group), colnames(norm_matrix)))
  norm_matrix <- norm_matrix[,shared_id]
  norm_matrix <- norm_matrix[,match(shared_id, colnames(norm_matrix))]
  sample_group <- sample_group[match(shared_id, rownames(sample_group)),]
  
  if ( is.null(gene_file) ) {
    return(list('count'=norm_matrix, 'group'=sample_group))
  } else {
    gene_list <- read.delim(gene_file, header=F, stringsAsFactors=F, check.names=F)
    return(list('count'=norm_matrix, 'group'=sample_group, 'gene'=gene_list))
  }
}



BoxPlotForGene <- function ( matrix, sample_group, gene_list ) {
  for ( gene in gene_list ) {
    if ( gene %in% rownames(matrix) ) {
      df <- data.frame(t(norm_matrix[gene,]), sample_group$Group)
      colnames(df) <- c("Counts", "Group")
      print(
        ggboxplot(df, x="Group", y="Counts", add="jitter") + 
          stat_compare_means() + 
          labs(title=gene, x="Group", y="Normalized Count") +
          theme(text=element_text(size=10), axis.text=element_text(size=10), 
                axis.title=element_text(size=12), plot.title=element_text(size=16), 
                legend.title=element_text(size=12), legend.text=element_text(size=12), 
                plot.subtitle=element_text(size=10), plot.caption=element_text(size=10))
      )
    } else {
      print(paste(gene, "not in matrix"))
    }
    readline(prompt="Press [enter] to continue")
  }
}



read_mat <- ReadFiles ( matrix_file=matrix_file, sample_group_file=sample_group_file, gene_file=gene_file)
norm_matrix <- read_mat$count
sample_group <- read_mat$group
ref_genes <- read_mat$gene

BoxPlotForGene ( norm_matrix, sample_group, ref_genes$V1 )



