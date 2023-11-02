#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @ChangeLog
#     20230222, only remove batch effect and normalization


suppressMessages(library(optparse))
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(limma))
suppressMessages(library(pheatmap))
suppressMessages(library(stringr))


option_list = list(
  make_option("--exp", type="character", default=NULL, help="expression matrix", metavar="character"),
  make_option("--group", type="character", default=NULL, help="sample group file", metavar="character"),
  make_option("--out", type="character", default=NULL, help="output prefix with directory", metavar="character"),
  make_option("--refgenes", type="character", default=NULL, help="reference gene list to normalize counts, either gene names separated by comma',' OR a file with each gene in one row", metavar="character"),
  make_option("--fit", type="character", default=NULL, help="saved lm fit model for cofactors", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



categorical_factor <- c('Group')
numerical_factor <- c('Age', 'RIN', 'Lib_conc', 'Data_amount')



RefGeneParse <- function ( refgenes ) {
  if ( file.exists(refgenes) ) {
    refgenes <- read.table(refgenes)$V1
  } else {
    # refgenes <- "PF4,PPBP,GP5,ITGA2B,TUBB1,SPARC,MYL9,GNG11"
    refgenes <- unlist(strsplit(refgenes, split=","))
  }
  return(refgenes)
}



NormalizeCounts <- function( df, refgenes=NULL ) {
  if ( is.null(refgenes) ) {
    # https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
    normalize_factor <- estimateSizeFactorsForMatrix(df)
  } else {
    # https://support.bioconductor.org/p/115682/
    # https://www.biostars.org/p/9522393/
    refgenes <- RefGeneParse(refgenes)
    refgenes <- refgenes[refgenes %in% rownames(df)]
    print(paste("Ref genes found:", paste(refgenes, collapse = ' ')))
    print(paste("Ref genes NOT found:", paste(refgenes[! refgenes %in% rownames(df)], collapse = ' ')))
    refgenes_idx <- which(rownames(df) %in% refgenes) #match(refgenes, rownames(df))
    normalize_factor <- estimateSizeFactorsForMatrix(df, controlGenes=refgenes_idx)
    if ( all(is.na(normalize_factor)) ) {
      normalize_factor <- estimateSizeFactorsForMatrix(df)
    }
  }
  normalized_counts <- df / normalize_factor
  normalized_counts <- log2(normalized_counts+1)
  
  return(normalized_counts)
}



removeBatchEffectWithFit <- function (x, batch = NULL, batch2 = NULL, covariates = NULL, 
                               design = matrix(1, ncol(x), 1), fit=NULL) 
{
  if (is.null(batch) && is.null(batch2) && is.null(covariates)) 
    return(as.matrix(x))
  if (!is.null(batch)) {
    batch <- as.factor(batch)
    contrasts(batch) <- contr.sum(levels(batch))
    batch <- model.matrix(~batch)[, -1, drop = FALSE]
  }
  if (!is.null(batch2)) {
    batch2 <- as.factor(batch2)
    contrasts(batch2) <- contr.sum(levels(batch2))
    batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
  }
  if (!is.null(covariates)) 
    covariates <- as.matrix(covariates)
  X.batch <- cbind(batch, batch2, covariates)
  if ( is.null(fit) ) {
    fit <- lmFit(x, cbind(design, X.batch))
  } else {
    load(fit)
  }
  beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
  beta[is.na(beta)] <- 0
  new_matrix <- as.matrix(x) - beta %*% t(X.batch)
  return(list('matrix'=new_matrix, 'fit'=fit))
}



RemoveCofactorEffect <- function ( df, sample_group, group='Group', fit=NULL ) {
  
  # removeBatchEffect: This simply removes any shifts in the log2-scale expression data that can be explained by batch. The design argument is necessary to avoiding removing variation associated with the treatment conditions.
  
  mm <- model.matrix(as.formula(paste0('~',group)), sample_group)
  factors <- colnames(sample_group)
  print(factors)
  batch_ids <- setdiff(factors[factors %in% categorical_factor], group)
  if ( length(batch_ids) == 0 ) {
    batch <- NULL
  } else {
    batch <- as.data.frame(sample_group[,batch_ids])
    colnames(batch) <- batch_ids
    for ( i in batch_ids ) {
      unique_level <- unique(batch[,i])
      replace_mat <- cbind(unique_level, index=range(0:(length(unique_level)-1)))
      for ( j in 1:nrow(replace_mat) ) {
        batch[,i] <- str_replace(batch[,i], replace_mat[j, 'unique_level'], replace_mat[j, 'index'])
      }
    }
  }
  covariates_list <- setdiff(numerical_factor, group)
  rmBatchlist <- removeBatchEffectWithFit(df,
                               batch = batch,
                               covariates = sample_group[, factors[factors %in% covariates_list]],
                               design = mm, fit=fit)

  return(list('matrix'=rmBatchlist$matrix, 'fit'=rmBatchlist$fit))
  
}



DrawHeatmap <- function (df, annoDF, out) {
  
  annotation_col <- as.data.frame(annoDF)
  
  pdf(out, width = ncol(df)/3, height = 10)
  print(pheatmap(df, #scale="row", 
                 cluster_rows=T, show_rownames=F,
                 cluster_cols=T, show_colnames=T, annotation_col=annotation_col))
  dev.off()
  
}



Main <- function ( matrix_file, sample_group_file, out, refgenes=NULL, group='Group', fit=NULL) {
  
  # read count file
  df <- read.delim(matrix_file, header=T, row.names=1, sep="\t", stringsAsFactors=F, check.names=F)
  sample_group <- read.delim(sample_group_file, header=T, row.names=1, sep="\t", stringsAsFactors=F, check.names=F)
  
  # change NAs to 0 and filter matrix
  sample_group_withNA <- sample_group
  sample_group[is.na(sample_group)] <- 0
  
  # shared samples
  shared_id <- sort(intersect(rownames(sample_group), colnames(df)))
  df <- df[,shared_id]
  df <- df[,match(shared_id, colnames(df))]
  sample_group <- sample_group[match(shared_id, rownames(sample_group)),]
  

  outdir <- dirname(out)
  if ( ! dir.exists(outdir) ) {
    dir.create(file.path(outdir), recursive = TRUE)
  }
  
  # normalization
  normalized_counts <- NormalizeCounts( df, refgenes=refgenes )
  
  # remove cofactor effect
  rmBatch <- RemoveCofactorEffect(normalized_counts, sample_group=sample_group, group=group, fit=fit)
  rmBatch$matrix[ rmBatch$matrix < 0 ] <- 0
  write.table(rmBatch$matrix, file=paste(out,"RemoveCofactorEffect","matrix","tsv",sep="."), append=F, quote=F, row.names=T, col.names=NA, sep="\t")
  if ( is.null(fit) ) {
    fit <- rmBatch$fit
    save(fit, file = paste(out,"RemoveCofactorEffect","fit","RData",sep=".") )
  }
  
  # heatmap, only use genes with all sample have > 0 value
  remove_zero_exp <- rmBatch$matrix[rowSums(rmBatch$matrix==0) == 0,]
  DrawHeatmap(df=remove_zero_exp, 
              annoDF=sample_group_withNA[, colnames(sample_group_withNA) %in% c(categorical_factor, numerical_factor)], 
              out=paste(out,"RemoveCofactorEffect","heatmap","pdf",sep=".")) 
  
}


Main( matrix_fil=opt$exp, sample_group_file=opt$group, out=opt$out, refgenes=opt$refgenes, group='Group', fit=opt$fit) 
# Main( matrix_file = 'merged.htseq.count.merge_two_batch', sample_group_file='conditions.with_cofactor.all_sample_two_batch.txt', out='test', refgenes=NULL, group='Group', fit=NULL)
# Main( matrix_file = 'merged.htseq.count', sample_group_file='condition.with_cofactor.all_samples.txt', out='test_deseq//test.all_sample', refgenes=NULL, group='Group', fit="test_deseq//test.RemoveCofactorEffect.fit.RData") 
