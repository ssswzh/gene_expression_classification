#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @ChangeLog
#     20230214, first version, adapted from zhangshouwei conduct_DESeq2_analysis.R
#     20230221, add removeBatchEffect DrawHeatmapByFeature


suppressMessages(library(optparse))
suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
# suppressMessages(library(ggbiplot))
suppressMessages(library(limma))
suppressMessages(library(gplots))
suppressMessages(library(pheatmap))


usage <- "
Usage:
    Rscript --vanilla conduct_DESeq2_analysis.modified.R --exp merged.count --group sample_group --compare leiomyoma,sarcoma --out out --dir ./

compare: first group is control, second group is test

"


# define arguments
option_list = list(
  make_option("--exp", type="character", default=NULL, help="expression matrix", metavar="character"),
  make_option("--group", type="character", default=NULL, help="sample group file", metavar="character"),
  make_option("--compare", type="character", default=NULL, help="compare groups, seperate by comma", metavar="character"),
  make_option("--out", type="character", default=NULL, help="output prefix", metavar="character"),
  make_option("--dir", type="character", default='./', help="output outdir, default current path", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



if ( FALSE ) {
  setwd("F://项目//TEP子宫肉瘤//tpm_matrix//merge_batch1_batch2//")
  exp_mat_file <- "merged.htseq.count.merge_two_batch" ## merged.htseq.count
  sample_group_file <- "conditions.with_cofactor.all_sample_two_batch.txt" ##condition.txt, which include sample information
  compare <- "Leiomyoma,Sarcoma"
  cmp = unlist(strsplit(compare,',')) ## compare.txt, wich include two subtypes, we put control on subtype1
  print(cmp)
  outdir <- 'test_deseq'
  dir.create(file.path(outdir), recursive = TRUE, showWarnings = F)
}


cofactor_list <- list(c('Batch','Group'),
                      c('Age','Batch','Group'),
                      c('Age','Lib_conc','Extraction_batch','Batch','Group'),
                      c('Age','Lib_conc','Extraction_batch','Data_amount','Batch','Group'),
                      c('Age','Hemolysis','RIN','Lib_conc','Extraction_batch','Data_amount','Batch','Group'))


categorical_factor <- c('Hemolysis', 'Extraction_batch', 'Batch', 'Group')
numerical_factor <- c('Age', 'RIN', 'Lib_conc', 'Data_amount')


ColumnAsFactor <- function ( sample_group ) {
  
  sample_group_factor <- sample_group
  for ( factor in colnames(sample_group) ) {
    if ( factor == 'Age' ) {
      sample_group_factor[, factor] <- cut(sample_group_factor[, factor], breaks = c(0,20,40,60,80,100))
      sample_group_factor[, factor] <- relevel(sample_group_factor[, factor], "(0,20]")
    }
    if ( factor == 'RIN' ) {
      sample_group_factor[, factor] <- cut(sample_group_factor[, factor], breaks = c(0,3,7,Inf))
      sample_group_factor[, factor] <- relevel(sample_group_factor[, factor], "(0,3]")
    }
    if ( factor == 'Lib_conc' ) {
      sample_group_factor[, factor] <- cut(sample_group_factor[, factor], breaks = c(0,20,40,60,80,Inf))
      sample_group_factor[, factor] <- relevel(sample_group_factor[, factor], "(0,20]")
    }
    if ( factor == 'Data_amount' ) {
      sample_group_factor[, factor] <- cut(sample_group_factor[, factor], breaks = c(0,10,20,Inf))
      sample_group_factor[, factor] <- relevel(sample_group_factor[, factor], "(0,10]")
    }
  }
  
  return(sample_group_factor)
  
}



PreprocessMatrix <- function ( exp_mat, sample_group ) {

  # filter genes
  # exp_mat$median <- apply(exp_mat[,-1], 1, median)
  # exp_mat <- exp_mat[order(exp_mat$name, exp_mat$median, decreasing=T),]
  # exp_mat <- exp_mat[,-c(1,ncol(exp_mat))]
  exp_mat <- exp_mat[!duplicated(exp_mat$name),]
  rownames(exp_mat) <- exp_mat$name
  exp_mat <- exp_mat[,-c(1)]
  
  # exclude genes less than five read counts in all samples, exclude samples with gene number < 5000
  sample_gene3k <- colnames(exp_mat)[colSums(exp_mat > 0) > 5000]
  exp_mat <- exp_mat[,sample_gene3k]
  exp_mat <- exp_mat[rowSums(exp_mat<5)!=ncol(exp_mat),]
  exp_mat <- exp_mat[rowSums(exp_mat)>=30,]
  exp_mat <- exp_mat[apply(exp_mat, 1, sd)>=1,]
  
  # only include samples in sample file and reorder
  shared_id <- sort(intersect(rownames(sample_group), colnames(exp_mat)))
  exp_mat <- exp_mat[,shared_id]
  exp_mat <- exp_mat[,match(shared_id, colnames(exp_mat))]
  sample_group <- sample_group[match(shared_id, rownames(sample_group)),]
  
  return(list('matrix'=exp_mat, 'group'=sample_group))
  
}



# DrawHeatmapByFeature <- function (DESeq2_object, out, genes) {
#   
#   # select <- order(rowMeans(counts(DESeq2_object,normalized=TRUE)),decreasing=TRUE)[1:genes]
#   annotation_col <- as.data.frame(colData(DESeq2_object)[,-c(ncol(colData(DESeq2_object)))])
#   
#   pdf(out, width = ncol(DESeq2_object)/3, height = max(length(genes)/2,10))
#   print(pheatmap(assay(normTransform(DESeq2_object))[genes,], cluster_rows=T, show_rownames=T,
#            cluster_cols=T, annotation_col=annotation_col))
#   dev.off()
#   
# }



BuildTestDEseq2Genes <- function (exp_mat, sample_group, formula_cols, contrast, out, reduced_cols=NULL, group='Group') {
  
  sample_group <- sample_group[,formula_cols]
  nonNAsamples <- rownames(sample_group[rowSums(is.na(sample_group))==0,])
  
  cofactor <- paste(formula_cols, collapse='+')
  formula <- paste('~', cofactor)
  sample_group_factor <- ColumnAsFactor(sample_group)
  DESeq2_data <- DESeqDataSetFromMatrix(countData = exp_mat[,nonNAsamples],
                                        colData = sample_group_factor[nonNAsamples,],
                                        design = as.formula(formula))
  
  # Extract results of DESeq2 analysis
  
  if ( is.null(reduced_cols) ) {
    DESeq2_object <- DESeq(DESeq2_data, quiet=T)
  } else {
    # http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test
    reduced <- paste('~', paste(reduced_cols, collapse='+'))
    DESeq2_object <- DESeq(DESeq2_data, test="LRT", reduced=as.formula(reduced), quiet=T)
  }
  
  DESeq2_result <- results(DESeq2_object, contrast=contrast)
  DESeq2_result <- DESeq2_result[order(DESeq2_result$padj),]
  DESeq2_result_Sig <- subset(DESeq2_result, padj<0.05)
  
  outfile <- paste(out, cofactor, sep=".")
  print(outfile)
  print(dim(sample_group))
  summary(DESeq2_result)
  write.table(as.matrix(DESeq2_result), file=paste(outfile,"DEseqResult","tsv",sep="."), append=F, quote=F, row.names=T, col.names=NA, sep="\t")
  
  # # get top genes for heatmap
  # if ( nrow(DESeq2_result_Sig) >10 ) {
  #   topgenes <- rownames(DESeq2_result_Sig)[1:20]
  # } else {
  #   new_set <- subset(DESeq2_result, padj<0.1)
  #   new_set <- new_set[order(new_set$padj),]
  #   topgenes <- rownames(new_set)[1:20]
  # }
  # topgenes <- topgenes[!is.na(topgenes)]
  # if ( length(topgenes) >= 2 ) {
  #   DrawHeatmapByFeature(DESeq2_object=DESeq2_object, out=paste(outfile,paste0('top',length(topgenes)),"heatmap","pdf",sep="."), genes=topgenes)
  # }
  
  return(list('object'=DESeq2_object, 'result'=DESeq2_result))
  
}



Main <- function ( exp_mat_file, sample_group_file, out, compare=NULL, outdir='./' ) {
  
  cmp = unlist(strsplit(compare,',')) ## compare.txt, wich include two subtypes, we put control on subtype1
  print(cmp)
  if ( ! dir.exists(outdir) ) {
    dir.create(file.path(outdir), recursive = TRUE)
  }
  out <- paste(outdir, out, sep='/')
  
  # Load data
  exp_mat <- read.delim(exp_mat_file, header=T, sep="\t", stringsAsFactors=F,check.names=F)
  sample_group <- read.delim(sample_group_file,header=T,row.names=1,sep="\t", stringsAsFactors=F,check.names=F)
  
  ## Preprocess data
  ## filter groups
  if ( 'Group' %in% colnames(sample_group) ) {
    group <- 'Group'
  } else if ( 'subtype' %in% colnames(sample_group) ) {
    group <- 'subtype'
  }
  
  sample_group <- sample_group[sample_group[,group] %in% cmp,]
  for ( fac in categorical_factor ) { 
    if ( fac %in% colnames(sample_group) ) {
    sample_group[,fac] <- factor(sample_group[,fac], levels = unique(sample_group[,fac]) )
    # sample_group[,fac] <- relevel(sample_group[,fac], ref = unique(sample_group[,fac])[1])
    }
  }
  
  preprocessed_result <- PreprocessMatrix(exp_mat, sample_group)
  filtered_mat <- preprocessed_result$matrix
  sample_group <- preprocessed_result$group
  
  for ( i in cofactor_list ) {
    reduced_cols <- setdiff(i, group)
    de <- BuildTestDEseq2Genes(filtered_mat, sample_group, formula_cols=i, contrast=c(group, cmp[2], cmp[1]), out=out, reduced_cols=reduced_cols)
    
  }
  
  
}


Main(exp_mat_file=opt$exp, sample_group_file=opt$group, out=opt$out, compare=opt$compare, outdir=opt$dir)
# Main(exp_mat_file="merged.htseq.count.merge_two_batch", sample_group_file="conditions.with_cofactor.all_sample_two_batch.txt", out="test", compare="Leiomyoma,Sarcoma", outdir="./")



# for ( i in 1:ncol(sample_group) ) {
#   sample_group[,i] <- as.character(sample_group[,i])
# }
# 
# design_model <- model.matrix(~ -1 + Age+Hemolysis+RIN+Lib_conc+Extraction_batch+Data_amount+Group, sample_group)
# 
# for ( i in 1:ncol(design_model) ) {
#   type <- unique(design_model[,i])
#   print(type)
#   print(design_model[,i])
#   design_model[,i] <- factor(design_model[,i], levels=type)
#   print(design_model[,i])
# }


# DESeq2_data <- DESeqDataSetFromMatrix(countData = filtered_mat[,rownames(design_model)],
#                                       colData = sample_group,
#                                       design = ~ Group)
# 
# DESeq2_object <- DESeq(DESeq2_data, quiet=T)


# library("AnnotationDbi")
# library("org.Hs.eg.db")
# 
# DESeq2_result$ENSEMBL <- mapIds(org.Hs.eg.db,
#                      keys=row.names(DESeq2_result),
#                      column="ENSEMBL",
#                      keytype="SYMBOL",
#                      multiVals="first")
# DESeq2_result$entrez <- mapIds(org.Hs.eg.db,
#                      keys=row.names(DESeq2_result),
#                      column="ENTREZID",
#                      keytype="SYMBOL",
#                      multiVals="first")

