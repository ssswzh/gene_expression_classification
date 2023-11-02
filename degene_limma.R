#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @Reference: 
#     https:\\bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html
#     https:\\www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
#     https:\\github.com/xuzhougeng/Learn-Bioinformatics/blob/master/7.Differential-Expression-Analysis.md
#     https:\\ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
#     https:\\www.jieandze1314.com/post/cnposts/90/
# @ChangeLog
#     20230227, first attempt
#     20230301, finished first version, ReadFiles(), SampleFilter(), LogCPMPlot(), NormFactorsPlot(), MDSforMatrix(), PCAforMatrix()
#     20230302, add PerformDEanalysis(), FillNAvalues()
#     20230418, add GeneFilter(), SampleCorrelation(), tSNEforMatrix(), BoxplotArrange()
#     20230419, add CheckConditions(), Gridggarrange(), BuildDesignModel(), DEGbyCofactor(), change BoxplotArrange() to BoxplotforMatrix()
#     20230420, add save=FALSE for functions for plotting
#     20230421, add OutputPrefix(), FeatureCorrelation(), Preprocess(), Normalization(), ClinicalTesting(), Correction(), Main() and options



suppressMessages(library(optparse))
suppressMessages(library(limma))
library(edgeR)
library(stringr)



if ( FALSE ) {
  setwd("F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\")
  exp_mat_file <- "merged.htseq.raw.count.merge_two_batch_used"
  sample_group_file <- "conditions.with_cofactor.all_sample_two_batch.addHealthy.tsv" 
  gene_len_file <- NULL
  outdir <- OutputPrefix( rootdir="limma_RNA_healthy", outdir=NULL, outprefix=NULL )
}



OutputPrefix <- function ( rootdir=NULL, outdir=NULL, outprefix='output' ) {
  if ( is.null(rootdir) ) {
    rootdir <- getwd()
  }
  if ( is.null(outdir) ) {
    outdir <- ""
  }
  outdir <- paste0(gsub("\\$", "", rootdir), '\\', outdir)
  if( !dir.exists(outdir) ) {
    print(paste("Create directory:", outdir))
    dir.create(outdir, recursive=T)
  }
  if ( ! is.null(outprefix) ) {
    outprefix <- paste0(gsub("\\$", "", outdir), '\\', outprefix)
    return(outprefix)
  } else {
    return(outdir)
  }
}



ReadFiles <- function ( exp_mat_file, sample_group_file, gene_len_file=NULL ) {
  
  exp_mat <- read.delim(exp_mat_file, header=T, sep="\t", row.names=1, stringsAsFactors=F,check.names=F)
  print( paste0("Raw matrix dimension: ", "(", nrow(exp_mat), ",", ncol(exp_mat), ")") )
  
  sample_group <- read.delim(sample_group_file, header=T, row.names=1, sep="\t", stringsAsFactors=F, check.names=F)
  print( paste0("Raw sample group dimension: ", "(", nrow(sample_group), ",", ncol(sample_group), ")") )
  
  # only include samples in sample file and reorder
  shared_id <- sort(intersect(rownames(sample_group), colnames(exp_mat)))
  exp_mat <- exp_mat[,shared_id]
  exp_mat <- exp_mat[,match(shared_id, colnames(exp_mat))]
  sample_group <- sample_group[match(shared_id, rownames(sample_group)),]
  print( paste0("Intersection of matrix and group: ", "(", nrow(exp_mat), ",", ncol(exp_mat), ")") )
  
  if ( is.null(gene_len_file) ) {
    return(list('count'=exp_mat, 'group'=sample_group))
  } else {
    gene_len <- read.delim(gene_len_file, header=T, sep="\t", stringsAsFactors=F, check.names=F)
    return(list('count'=exp_mat, 'group'=sample_group, 'gene'=gene_len))
  }
  
}



FillNAvalues <- function ( df, by='Group', fill='mean' ) {
  
  df_order <- rownames(df)
  new_df <- data.frame(matrix(nrow=0,ncol=ncol(df)))
  colnames(new_df) <- colnames(df)
  groups <- unique(df[,by])
  cols2fill <- setdiff(colnames(df), by)
  for ( grp in groups ) {
    tmpdf <- df[df[,by]==grp,]
    for ( col in cols2fill ) {
      if ( fill == 'median' ) { fill_value <- median(tmpdf[,col], na.rm=TRUE) }
      if ( fill == 'mean' ) { fill_value <- mean(tmpdf[,col], na.rm=TRUE) }
      tmpdf[is.na(tmpdf[,col]), col] <- fill_value
    }
    new_df <- rbind(new_df, tmpdf)
  }
  new_df <- new_df[match(df_order, rownames(new_df)),]
  return(new_df)
  
}



GeneExpressionPlot <- function ( counts, highlight_genes=NULL, cal='median', save=FALSE ) {
  # overall gene expression levels across all samples
  # highlight_genes: number or vector
  # cal: 'sum', 'mean', 'median'
  
  if ( ! cal %in% c('sum', 'mean', 'median') ) {
    stop("'cal' should be one of c('sum', 'mean', 'median')")
  }
  if ( (! is.null(highlight_genes)) && (! typeof(highlight_genes)=="character") && (! typeof(highlight_genes)=="double") ) {
    stop("'highlight_genes' should be a number or a vector")
  }
  
  library(scales)
  expDf <- data.frame('GeneExp'=sort(apply(counts, 1, median), decreasing = T), 'Rank'=seq(1,nrow(counts)))
  expDf$Gene <- rownames(expDf)
  
  # basic plot and if log coordinate
  plt <- ggplot(expDf, aes(x=Rank, y=GeneExp)) + geom_point(col='darkgreen') + theme_classic()
  if ( max(expDf$GeneExp) > 1000 ) {
    plt <- plt + scale_y_continuous(trans = log10_trans(),
                                    breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x)))
  }
  
  # highlight genes
  if ( ! is.null(highlight_genes) ) {
    if ( typeof(highlight_genes)=="character" ) {
      plt <- plt + geom_text(aes(label=ifelse(Gene %in% highlight_genes, as.character(Gene), '')), col='red', hjust=0, vjust=0)
    } else if ( typeof(highlight_genes)=="double" ) {
      plt <- plt + geom_text(aes(label=ifelse(GeneExp>highlight_genes, as.character(Gene),'')), col='red', hjust=0, vjust=0)
    } 
  }
  
  if ( save != FALSE ) {
    pdf(paste(save,"GeneExp",cal,'pdf',sep="."), width=6, height=6)
    print(plt)
    dev.off()
  } else {
    print(plt)
  }
  
  return(expDf)
  
}



GeneFilter <- function ( counts, number=10, ratio=0.1 ) {
  
  # exclude genes with <number reads in at least samples*ratio
  genes <- rownames(counts)[rowSums(counts < number) > ncol(counts)*ratio]
  counts <- counts[! rownames(counts) %in% genes, ]
  print( paste("Left genes:", nrow(counts)) )
  return(list('counts'=counts, 'genes'=genes))
  
}



SampleFilter <- function ( counts, sample_group, gene_count=5, gene_num=5000 ) {
  
  # 样本质控，表达（count>5）的基因数量>1500的样本
  keep_samples <- colnames(counts)[colSums(counts > gene_count) > gene_num]
  discard <- setdiff(colnames(counts), keep_samples)
  sample_group <- sample_group[keep_samples,]
  counts <- counts[,keep_samples]
  print(paste("Left samples:", ncol(counts)) )
  print(paste("Discard samples:", length(discard)) )
  print(paste(discard, collapse=" "))
  return(list('counts'=counts, 'group'=sample_group, 'discard'=discard))
  
}



SampleCorrelation <- function ( counts, sample_group=NULL, by=NULL, cutoff=0.4, cal='sum', method='pearson', save=FALSE ) {
  # 根据样本相关性过滤样本, leave one out
  # cal = c('sum', 'mean', 'median', NULL)
  # method = c("pearson", "kendall", "spearman")
  # save = FALSE or path to save figure
  
  if ( is.null(sample_group) ) {
    print("SampleCorrelation: no sample info is provided, ignore option 'by'")
  }
  if ( (! is.null(by)) && (! by %in% colnames(sample_group)) ) {
    stop(paste0('ERROR: by=',by,'not in sample_group'))
  }
  if ( ! (cal %in% c('sum', 'mean', 'median') || is.null(cal)) ) {
    stop(paste0('ERROR: cal=',cal,'not in supported, \nvalid choices are: ', paste(c('sum', 'mean', 'median', 'NULL'), collapse=" ")))
  }
  
  # make up a fake group
  if ( is.null(by) ) {
    by <- 'All'
    sample_group[,by] <- 'All'
  }
  
  # calculate correlation by each group
  dfList <- list()
  correlationDf <- data.frame()
  groups <- unique(sort(sample_group[,by]))
  for (group in groups) {
    groupsamples <- rownames(sample_group[sample_group[,by] %in% group, ])
    groupDf <- counts[,groupsamples]
    dfList[[group]] <- data.frame(cor(groupDf, method=method), check.names = F)
    if ( ! is.null(cal) ) { # if call==NULL, do not perform leave-one-out correlation
      for ( s in colnames(groupDf) ) {
        sampleDf <- groupDf[,s]
        tmpDf <- groupDf[,! colnames(groupDf) %in% s]
        tmpDf <- apply(tmpDf, 1, cal)
        score <- cor(tmpDf, sampleDf, method=method)
        correlationDf[s,'Score'] <- score
        correlationDf[s,'Group'] <- group
      }
    }
  }
  
  # plot
  if ( save != FALSE ) {
    # heatmap
    save <- paste(save, 'sample_correlation', by, method, paste0(cal,cutoff), sep='.')
    library(pheatmap)
    for ( ll in names(dfList) ) {
      mat <- dfList[[ll]]
      pheatmap(mat, display_numbers=T, cluster_rows=F, cluster_cols=F,
               width=max(12, ncol(mat)/4), height=max(12, ncol(mat)/4), 
               filename=paste(save, 'correlation_heatmap', ll, 'pdf', sep='.'))
      mat$Mean <- (rowSums(mat)-1)/(ncol(mat)-1)
      write.table(mat, paste(save, 'correlation_heatmap', ll, 'tsv', sep='.'), sep="\t", col.names=NA, quote=F)
    }
    # boxplot if cal!=NULL
    if ( ! is.null(cal) ) {
      library(ggplot2)
      library(ggpubr)
      lv <- unique(sample_group[,by])
      sample_group[,by] <- factor(sample_group[,by], levels=lv)
      pdf(paste(save, 'correlation_score', 'pdf', sep='.'), width=max(6, length(unique(correlationDf$Group))), height=8)
      if ( length(lv) == 1 ) {
        print(ggplot(correlationDf, aes(x=Group, y=Score, color=Group)) + 
                geom_hline(yintercept=cutoff, color='red', linetype='dashed') + 
                geom_boxplot(outlier.color=NA) + geom_point(position="jitter") + theme_classic())
      } else {
        comparisons <- as.list(data.frame(combn(unique(as.character(sample_group[,by])), 2)))
        print(ggplot(correlationDf, aes(x=Group, y=Score, color=Group)) + 
                stat_compare_means(comparisons = comparisons) + 
                geom_hline(yintercept=cutoff, color='red', linetype='dashed') + 
                geom_boxplot(outlier.color=NA) + geom_point(position="jitter") + theme_classic())
      }
      dev.off()
      write.table(correlationDf, paste(save,'correlation_score','tsv',sep='.'), sep='\t', col.names=NA, quote=F)
    }
  }
  
  # exclude samples in matrix
  exclude_samples <- correlationDf[correlationDf$Score < cutoff, ]
  if ( nrow(exclude_samples) > 20 ) {
    print( paste("Exclude samples:", nrow(exclude_samples)) )
    print(paste(rownames(exclude_samples), collapse=" "))
  } else {
    print("Exclude samples:")
    print(exclude_samples)
  }
  counts <- counts[, ! colnames(counts) %in% exclude_samples]
  sample_group <- sample_group[! sample_group %in% exclude_samples,]
  
  return_list <- list('count'=counts, 'group'=sample_group, 'correlation'=correlationDf, 'discard'=exclude_samples)
  for ( i in names(dfList) ) { return_list[[i]] <- dfList[[i]]}
  return(return_list)
  
}



LogCPMPlot <- function ( dge_counts, dge_counts_filtered, sampleName, geneLength=NULL, save=FALSE ) {
  # count转化为log2的值，log-CPM值是根据CPM值通过log2(CPM + prior.count/L)计算得到的，L是样本文库大小（以百万计）的平均值
  
  if ( is.null(geneLength) ) {
    lcpm <- cpm(dge_counts, log=TRUE, prior.count=2)
    lcpm_filtered <- cpm(dge_counts_filtered, log=TRUE, prior.count=2)
  } else {
    lcpm <- cpm(dge_counts, log=TRUE, prior.count=2, gene.length=geneLength)
    lcpm_filtered <- cpm(dge_counts_filtered, log=TRUE, prior.count=2, gene.length=geneLength)
  }
  
  # overall cutoff
  L <- mean(dge_counts$samples$lib.size) * 1e-6
  M <- median(dge_counts$samples$lib.size) * 1e-6
  lcpm.cutoff <- log2(10/M + 2/L)
  
  if ( save != FALSE ) {
    library(RColorBrewer)
    nsamples <- ncol(dge_counts)
    col <- rainbow(nsamples*1.2)
    
    pdf(paste(save, 'logCPM_distribution', 'pdf', sep='.'), width=16, height=8)
    par(mfrow=c(1,2))

    # raw
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
    title(main="Raw data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", sampleName, text.col=col, bty="n", cex = 0.5)

    # filtered
    plot(density(lcpm_filtered[,1]), col=col[1], lwd=2, ylim=c(0,0.6), las=2, main="", xlab="")
    title(main="Filtered data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
      den <- density(lcpm_filtered[,i])
      lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", sampleName, text.col=col, bty="n", cex = 0.5)
    
    dev.off()
  }
  
  return(list('cutoff'=lcpm.cutoff, 'logCPM'=lcpm, 'logCPM_filtered'=lcpm_filtered))
  
}



NormFactorsPlot <- function ( logCPM_Unnormalized, logCPM_normalized, save=FALSE ) {
  
  # 查看标准化前后的数据差异
  col <- rainbow(ncol(logCPM_Unnormalized)*1.2)
  if ( save != FALSE ) {
    pdf(paste(save, 'lib_normalization', 'pdf', sep='.'), width=max(20,ncol(logCPM_Unnormalized)/10), height=8)
  }
  par(mfrow=c(1,2))
  boxplot(logCPM_Unnormalized, las=2, col=col, main="")
  title(main="Unnormalised data",ylab="Log-cpm")
  boxplot(logCPM_normalized, las=2, col=col, main="")
  title(main="Normalised data",ylab="Log-cpm")
  
  if ( save != FALSE ) {
    dev.off()
  }
}



CheckConditions <- function ( dataframe, conditions ) {
  # check conditions in dataframe columns
  new_condition <- conditions[conditions %in% colnames(dataframe)]
  if ( length(new_condition) != length(conditions) ) {
    print("Names NOT in dataframe:")
    print(conditions[! conditions %in% colnames(dataframe)])
  }
  if ( length(new_condition) != 0 ) {
    return(new_condition)
  } else {
    stop(paste("No names in dataframe:", paste(conditions, collapse=" ")))
  }
}



Gridggarrange <- function ( list_plots, conditions, save=FALSE ) {
  # output ggarrange figures
  
  if ( save == FALSE ) { # open null device
    save <- NULL
  }
  
  if (length(conditions)/3 >= 1 ) {
    pdf(save, width=15, height=ceiling(length(conditions)/3)*5)
    plt <- ggarrange(plotlist = list_plots, common.legend = FALSE, 
                     ncol = 3, nrow = ceiling(length(conditions)/3))#, 
                     # width = 18, height = ceiling(length(conditions)/3)*4 )
    print( plt )
    dev.off()
  } else {
    pdf(save, width=length(conditions)*5, height=6)
    plt <- ggarrange(plotlist = list_plots, common.legend = FALSE, 
                     ncol = length(conditions), nrow = 1)#,
                     # width = length(conditions)*6, height = 6)
    print( plt )
    dev.off()
  }
  
  if ( is.null(save) ) { print( plt ) }
  
}



MDSforMatrix <- function ( matrix, sample_group, conditions=c('Hemolysis', 'Extraction_batch', 'Batch', 'Group'), save=FALSE ) {
  
  conditions <- CheckConditions( sample_group, conditions )
  colors <- list()
  for ( con in conditions ) {
    sample_group[,con] <- as.character(sample_group[,con])
    class <- unique(sample_group[,con])
    exp_cols <- brewer.pal(length(class), "Set1")
    fac_col <- sample_group[,con]
    for ( i in c(1:length(class)) ) {
      fac_col[fac_col==class[i]] <- exp_cols[i]
    }
    colors[[con]] <- fac_col
  }
  
  # plot
  if (length(conditions)/3 >= 1 ) {
    n_row <- ceiling(length(conditions)/3)
    n_col <- 3
    par(mfrow = c(n_row, 3))
  } else {
    n_row <- 1 
    n_col <- length(conditions)
  }
  
  if ( save != FALSE ) { # open null device
    save <- paste(save, 'MDS', 'pdf', sep='.')
    pdf(save, width=n_col*4, height=n_row*4)
  }
  
  par(mfrow = c(n_row, n_col))
  for ( con in conditions ) {
    plotMDS(matrix, label=sample_group[,con], col=colors[[con]], var.explained=T)
    title(con)
  }
  
  if ( save != FALSE ) {
    dev.off()
  }
  
}



PCAforMatrix <- function ( matrix, sample_group, conditions=c('Extraction_batch', 'Batch', 'Group'), save=FALSE ) {
  
  conditions <- CheckConditions( sample_group, conditions )
  
  library(FactoMineR)
  library(factoextra)
  new_pca <- PCA(t(matrix), graph = FALSE)
  
  for ( con in conditions ) {
    sample_group[,con] <- as.character(sample_group[,con])
  }

  # plot
  Map(function(con) {
    fviz_pca_ind(new_pca, geom.ind = "point",
                 col.ind = sample_group[,con], addEllipses=TRUE, legend.title=con, title=con)
  }, conditions) -> list_plots
  
  if ( save != FALSE ) { # open null device
    save <- paste(save, 'PCA', 'pdf', sep='.')
  }
  Gridggarrange( list_plots, conditions, save=save )
  
}



tSNEforMatrix <- function ( matrix, sample_group, conditions=c('Extraction_batch', 'Batch', 'Group'), save=FALSE ) {
  
  conditions <- CheckConditions( sample_group, conditions )
  
  library(Rtsne)
  library(ggplot2)
  sample_all <- t(matrix)
  tsne_out <- Rtsne(sample_all, pca=FALSE, perplexity=10, theta=0.0)
  # 获取tSNE的坐标值
  str(tsne_out)
  # 其中在Y中存储了绘制图坐标
  tsnes <- tsne_out$Y
  colnames(tsnes) <- c("tSNE1", "tSNE2") #为坐标添加列名
  # 在此基础上添加颜色分组信息，首先还是将tsnes这个矩阵变成数据框，然后增加一列group信息，最后映射在geom_point中
  tsnes=as.data.frame(tsnes)
  for ( con in conditions ) {
    oldfactor <- con
    con <- str_replace(con, " ", "_")
    tsnes[,con] <- as.character(sample_group[,oldfactor])
  }
  
  # plot
  Map(function(con) {
    ggplot(tsnes, aes(x = tSNE1, y = tSNE2)) + geom_point(aes_string(col=con)) + theme_classic() +
      stat_ellipse(aes_string(color=con), type = "norm")
  }, conditions) -> list_plots
  
  if ( save != FALSE ) { # open null device
    save <- paste(save, 'tSNE', 'pdf', sep='.')
  }
  
  Gridggarrange( list_plots, conditions, save=save )
  
}



BoxplotforMatrix <- function ( sample_group, by='Group', conditions=c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'), save=FALSE ) {
  
  conditions <- CheckConditions( sample_group, conditions )
  
  library(ggplot2, quietly=TRUE)
  theme_set(theme_classic())
  library(ggpubr)
  
  lv <- unique(sample_group[,by])
  sample_group[,by] <- factor(sample_group[,by], levels=lv)
  comparisons <- as.list(unname(data.frame( combn(unique(as.character(sample_group[,by])), 2) )))
  # print(comparisons)
  
  Map(function(col) {
    ggboxplot(sample_group, x=by, y=col, add="jitter", color=by) + 
      stat_compare_means(comparisons = comparisons) + 
      labs(title=col, x=by, y=col) +
      theme(text=element_text(size=10), axis.text=element_text(size=10), 
            axis.title=element_text(size=12), plot.title=element_text(size=16), 
            legend.title=element_text(size=12), legend.text=element_text(size=12), 
            plot.subtitle=element_text(size=10), plot.caption=element_text(size=10))
  }, conditions) -> list_plots

  if ( save != FALSE ) { # open null device
    save <- paste(save, "compare",by,"pdf",sep=".")
  }
  Gridggarrange( list_plots, conditions, save=save )
  
}



BuildDesignModel <- function ( sample_group, cofactors=c('Group', 'Age', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Batch', 'Conservation_days', 'Lib_amount', 'Data_amount') ) {
  
  design_factor <- colnames(sample_group)[colnames(sample_group) %in% cofactors]
  design_factor <- design_factor[match(cofactors, design_factor)]
  for(i in 1:length(design_factor)) {
    if(str_count(design_factor[i],' ')!=0){design_factor[i] <- paste0("`",design_factor[i],"`") }
  }
  formula <- paste0("~ 0 + ", paste0(design_factor, collapse=" + "))
  design_model <- model.matrix(as.formula(formula), sample_group)
  colnames(design_model) <- gsub("Group", "", colnames(design_model))
  
  return(list('model'=design_model, 'formula'=formula))
  
}



PerformDEanalysis <- function ( EList, design_model, contr.matrix=NULL, coefficients=NULL, out="limma.output", pvalue=0.05, logFC=1, plot=T ) {
  
  # limma的线性建模使用lmFit和contrasts.fit函数进行，每个基因的表达值都会单独拟合一个模型。然后通过借用全体基因的信息来进行经验贝叶斯调整（empirical Bayes moderation），这样可以更精确地估算各基因的差异性
  
  if ( (is.null(contr.matrix) && is.null(coefficients)) || (!is.null(contr.matrix) && !is.null(coefficients)) ) {
    stop("Please give one of contr.matrix and coefficients, but not both.")
  }
  
  # lmFit computes coefficients, residual variances and standard errors.
  vfit <- lmFit(EList, design_model)
  # head(coef(vfit),5)
  
  # contrasts.fit converts the coefficients and standard errors to reflect the contrasts rather than the original design matrix, but does not compute t-statistics or p-values.
  if  ( ! is.null(contr.matrix) ) {
    vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
  }
  if ( ! is.null(coefficients) ) {
    vfit <- contrasts.fit(vfit, coefficients=coefficients)
  }
  
  # empirical Bayes moderation, eBayes computes t-statistics and p-values from the coefficients and standard errors.
  efit <- eBayes(vfit)
  
  if ( plot ) {
    plotSA(efit, main="Final model: Mean-variance trend")
  }
  
  # 所有差异基因
  print(summary(decideTests(efit, adjust.method="fdr", p.value=pvalue, lfc=logFC)))
  degene_result <- topTable(efit, n = Inf, adjust = "fdr")
  degene_final = na.omit(degene_result)
  out <- gsub(" ","",out)
  out <- gsub("[*]","&",out)
  write.table(degene_final, paste(out,"tsv",sep="."), col.names = NA, row.names=T, sep="\t", quote=F)  
  
  ## 筛选出显著差异基因
  All_diffSig <- degene_final[(degene_final$adj.P.Val < pvalue & (degene_final$logFC>logFC | degene_final$logFC < (-logFC))),]
  new_col <- c('change', colnames(All_diffSig))
  All_diffSig$change <- ifelse(All_diffSig$adj.P.Val>pvalue, "NS", ifelse(All_diffSig$logFC>logFC, "UP", "DOWN") )
  All_diffSig <- All_diffSig[,match(new_col,colnames(All_diffSig))]
  write.table(All_diffSig, paste(out,"sig.tsv",sep="."), col.names = NA, row.names=T, sep="\t", quote=F)  
  
  return(list('DEG'=All_diffSig, 'efit'=efit))
}



FeatureCorrelation <- function ( matrix, out='out' ) {
  library(pheatmap)
  feature_correlation <- cor(t(matrix))
  if ( nrow(feature_correlation) > 50 ) {
    correlation_cluster <- pheatmap(round(feature_correlation,2), border_color=NA, #display_numbers=T,
                                    width=nrow(matrix)/3, height=nrow(matrix)/3, 
                                    filename=paste(out,"sig.feature_correlation.pdf",sep='.') )
  } else {
    correlation_cluster <- pheatmap(round(feature_correlation,2), display_numbers=T,
                                            width=nrow(matrix)/3, height=nrow(matrix)/3, 
                                            filename=paste(out,"sig.feature_correlation.pdf",sep='.') )
  }
  
  feature_correlation <- feature_correlation[rownames(feature_correlation[correlation_cluster$tree_row[["order"]],]), ]
  feature_correlation <- feature_correlation[, rownames(feature_correlation[,correlation_cluster$tree_col[["order"]]])]
  write.table(feature_correlation, paste(out,"feature_correlation.tsv",sep='.'), quote=F, sep="\t", col.names=NA)
  feature_mean <- data.frame('mean'=apply(matrix[rownames(feature_correlation),], 1, mean))
  write.table(feature_mean, paste(out,"feature_meanExp.tsv",sep='.'), quote=F, sep="\t", col.names=NA)
  
  return(list('correlation'=feature_correlation, 'meanExp'=feature_mean))
}



#  module functions

Preprocess <- function ( exp_mat_file, sample_group_file, outdir='outdir', 
                         high_conf_count=30, gene_filter=c(10,0.9), sample_filter=c(5,5000), fillNA=FALSE) {
  
  # read files 
  read_mat <- ReadFiles ( exp_mat_file, sample_group_file, gene_len_file=NULL)
  counts <- read_mat$count
  counts <- counts[rowSums(counts)!=0,]
  sample_group <- read_mat$group

  # clinical tests
  clinical_prefix <- OutputPrefix( rootdir=outdir, outdir='clinical_test', outprefix='count' )
  
  # 高置信基因在两组的差异
  high_conf_gene <- data.frame('High_confidence_Genes'=colSums(counts>=high_conf_count))
  high_conf_gene[,c('Group','Clinical_group')] <- sample_group[match(rownames(sample_group), rownames(high_conf_gene)), c('Group','Clinical_group')]
  for ( i in c('Group','Clinical_group') ) {
    BoxplotforMatrix( high_conf_gene, by=i, conditions=c('High_confidence_Genes'), save=paste(clinical_prefix,"high_confidence_genes",sep=".") )
  }
  
  # 样本间表达中位值最高的基因
  highlight_genes <- c('MT-RNR2', 'MT-RNR1', 'TMSB4X', 'B2M', 'MT_ND1')
  median_exp <- GeneExpressionPlot( counts, highlight_genes=highlight_genes, cal='median', save=clinical_prefix )
  
  # 过滤基因
  if ( ! is.null(gene_filter) || gene_filter!=FALSE ) {
    geneflt <- GeneFilter(counts, number=gene_filter[1], ratio=gene_filter[2])
    excluded_genes <- geneflt$genes
    counts <- geneflt$counts
  }
  
  # 过滤样本
  if ( ! is.null(sample_filter) || sample_filter!=FALSE ) {
    samflt <- SampleFilter(counts, sample_group, gene_count=sample_filter[1], gene_num=sample_filter[2])
    counts <- samflt$counts
    sample_group <- samflt$group
  }
  
  # 对不同分组的数据，如果有NA值则赋值组内平均值
  if ( fillNA ) {
    sample_group <- FillNAvalues( sample_group, by='Group', fill='mean')
  }
  
  # raw count matrix转换DEGList
  dge_counts <- DGEList(counts, samples=sample_group)
  dge_counts$samples <- dge_counts$samples[,-1]
  
  # filterByExpr函数提供了自动过滤基因的方法，可保留尽可能多的有足够表达计数的基因
  keep.exprs <- filterByExpr(dge_counts, group=sample_group$Group)
  dge_counts_filtered <- dge_counts[keep.exprs,, keep.lib.sizes=FALSE]
  keep_genes <- rownames(dge_counts_filtered)
  print( paste("Kept gene number by filterByExpr():", length(keep_genes)) )
  
  return(list('DGEList'=dge_counts, 'DGEListFiltered'=dge_counts_filtered, 'Group'=sample_group))
  
}



Normalization <- function ( dgelist, dgelist_filter, sample_group, save='outdir', 
                            category_group=c('Hemolysis', 'Extraction_batch', 'Reverse_batch', 'Batch', 'Clinical_group', 'Group') ) {
  
  # 过滤前后的logCPM密度图
  cpm_list <- LogCPMPlot( dgelist, dgelist_filter, sampleName=rownames(sample_group), geneLength=NULL, save=save )
  # logCPMcutoff <- cpm_list$cutoff
  # logCPM <- cpm_list$logCPM
  logCPM_filtered <- cpm_list$logCPM_filtered
  
  # 计算标准化因子并可视化标准化前后的logCPM分布
  dgelist_filter <- calcNormFactors(dgelist_filter, method = "TMM")
  dgelist_filter_logCPM <- cpm(dgelist_filter, log=TRUE)
  
  NormFactorsPlot(logCPM_filtered, dgelist_filter_logCPM, save=save)
  MDSforMatrix(dgelist_filter, sample_group, conditions=category_group, save=save )
  PCAforMatrix(dgelist_filter_logCPM, sample_group, conditions=category_group, save=save )
  tSNEforMatrix(dgelist_filter_logCPM, sample_group, conditions=category_group, save=save )
  # library(GLimma)
  # glMDSPlot(dgelist_filter_logCPM, labels=paste(sample_group$Group, sample_group$Batch, sep="_"), 
  # groups=dge_counts_filtered$samples[,category_group], launch=TRUE)
  
  return(list('DGEListFiltered'=dgelist_filter, 'logCPM'=dgelist_filter_logCPM) )
  
}



ClinicalTesting <- function ( dgelist, sample_group, outdir='outdir', 
                              numeric_cols=c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'), 
                              category_cols=c("Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Sequencing_batch")) {
  
  # 查看条件在两组中有无差异，如果有差异就不校正
  test_prefix <- OutputPrefix( rootdir=outdir, outdir='clinical_test', outprefix='sample_info' )
  
  # 数值型数据
  sample_group$Lib_size <- dgelist$samples[match(rownames(sample_group),rownames(dgelist$samples)),'lib.size']
  if ( is.null(numeric_cols) ) {
    numeric_cols <- c('Lib_size')
  } else {
    numeric_cols <- c(numeric_cols, 'Lib_size')
  }
  cmpby <- c('Group', 'Clinical_group')
  cmpby <- cmpby[cmpby %in% colnames(sample_group)]
  for( cmp in cmpby ) {
    BoxplotforMatrix( sample_group, by=cmp, conditions=numeric_cols, 
                      save=paste(test_prefix,"clinical_numeric",sep=".") )
  }
  # 在分组中无显著差异的条件（校正）：Age, RIN_value, cDNA_conc, Lib_conc, Data_amount, Lib_size
  # 在分组中有显著差异的条件（不校正）：Conservation_days, Lib_amount
  
  # 分类型数据
  if ( ! is.null(category_cols) ) {
    test_result <- data.frame()
    for ( i in category_cols ) {
      test <- fisher.test(table(sample_group[, c("Group", i)]), simulate.p.value=TRUE)
      test_result[i,"Pvalue"] <- test$p.value
    }
    write.table(test_result, paste(test_prefix, "clinical_category.fishertest.tsv",sep="."), col.names=NA, sep="\t", quote=F)
    for ( i in rownames(test_result)[test_result$Pvalue<=0.05] ) {
      print(table(sample_group[, c("Group", i)]))
    }
  }
  
}



Correction <- function ( dgelist, sample_group, outdir='outdir',
                         correction_factors=c('Group', 'cDNA_conc', 'Lib_conc', 'Data_amount', 'Batch', "Reverse_batch", "Sequencing_batch")) {
  
  # 在limma中，假设log-CPM值符合正态分布，因此在对RNA-seq的log-CPM值进行线性建模时，需要使用voom函数计算每个基因的权重从而调整均值与方差的关系，否则分析得到的结果可能是错误的。
  # voom()作用是原始counts转换为logCPM值，将所有计数加0.5，以避免取对数零。然后，将logCPM值矩阵进行标准化。
  # 根据每个基因的log2CPM制作了线性模型，并计算了残差，利用了平均表达量（红线）拟合了sqrt(residual standard deviation)；最后得到的平滑曲线可以用来得到每个基因和样本的权重
  # 如果横坐标接近0的位置出现迅速上升，说明low counts数比较多
  # 通常而言，方差是测序实验操作中的技术差异和来自不同细胞类群的重复样本之间的生物学差异的结合，而voom图会显示出一个在均值与方差之间递减的趋势。生物学差异高的实验通常会有更平坦的趋势，其方差值在高表达处稳定。生物学差异低的实验更倾向于急剧下降的趋势。
  
  voom_prefix <- OutputPrefix( rootdir=outdir, outdir='voom_correction', outprefix='voom_normalized' )
  
  # 构建校正数据的design model
  if ( 'Batch' %in% colnames(sample_group) ) {
    sample_group$Batch <- as.character(sample_group$Batch)
  }
  correction_model <- BuildDesignModel( sample_group, cofactors=correction_factors )
  design <- correction_model$model
  correction_formula <- correction_model$formula
  colnames(design) <- gsub("Group", "", colnames(design))
  print( paste0("Design model dimension: ", "(", nrow(design), ",", ncol(design), ")") )
  
  norm_exp <- voom(dgelist$counts[,rownames(design)], design, plot=F)
  write.table(norm_exp$E, paste(voom_prefix,"matrix",gsub(" ","",correction_formula),"tsv",sep="."), col.names = NA, row.names=T, sep="\t", quote=F)
  write.table(sample_group[colnames(norm_exp),], paste(voom_prefix,"sample_group",gsub(" ","",correction_formula),"tsv",sep="."), col.names = NA, row.names=T, sep="\t", quote=F)
  
  MDSforMatrix(norm_exp, sample_group[colnames(norm_exp),], conditions=correction_factors, save=voom_prefix )
  PCAforMatrix(norm_exp$E, sample_group[colnames(norm_exp),], conditions=correction_factors, save=voom_prefix )
  tSNEforMatrix(norm_exp$E, sample_group[colnames(norm_exp),], conditions=correction_factors, save=voom_prefix )

  # 未过滤
  # blank <- voom(dge_counts$counts[,rownames(design)], design, plot = T)
  
  return(norm_exp)
  
}



DEGbyCofactor <- function ( norm_exp, sample_group, cofactors, out='output' ) {
  # 构建检测差异基因的新design model
  build_model <- BuildDesignModel( sample_group, cofactors=cofactors )
  design_model <- build_model$model
  formula <- build_model$formula
  
  # 构建比对矩阵
  contr.matrix <- makeContrasts(Sarcoma - nonSarcoma, levels = design_model)
  
  group_degene_list <- PerformDEanalysis( EList=norm_exp[,colnames(norm_exp) %in% rownames(design_model)], 
                                          design_model=design_model, 
                                          contr.matrix=contr.matrix, 
                                          out=paste(out,formula,sep="."), 
                                          pvalue=0.05, logFC=1, plot=T )
  group_degene_list[['formula']] <- formula
  return(group_degene_list)
}




Main <- function ( count_file, group_file, 
                   outdir = 'outdir', 
                   gene_filter = FALSE, 
                   sample_filter = '5,5000', 
                   sample_cor = 'sum,pearson,0',
                   correction_factors = c('Group', "Age", "RIN_value"), 
                   de_factors = c('Group', "Age", "RIN_value", "Batch"),
                   numeric_cols = c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'),
                   category_cols = c("Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Sequencing_batch", 'Clinical_group', 'Group') ) {
  
  outdir <- OutputPrefix( rootdir=outdir, outdir=NULL, outprefix=NULL )
  
  # --- Read inputs and statistics --- #
  if ( ! (gene_filter==FALSE || is.null(gene_filter)) ) { gene_filter <- as.integer(unlist(str_split(gene_filter,','))) }
  if ( ! (sample_filter==FALSE || is.null(sample_filter)) ) { sample_filter <- as.integer(unlist(str_split(sample_filter,','))) }
  preprocess <- Preprocess( count_file, group_file, outdir=outdir, 
                            high_conf_count=30, gene_filter=gene_filter, sample_filter=sample_filter, fillNA=FALSE )
  sample_group <- preprocess$Group
  
  # --- Normalization and correlation --- #
  
  norm_prefix <- OutputPrefix( rootdir=outdir, outdir='normalization', outprefix='logCPM' )
  normalization <- Normalization( dgelist=preprocess$DGEList, dgelist_filter=preprocess$DGEListFiltered, sample_group=sample_group, 
                                  save=norm_prefix, category_group=category_cols )
  
  # 根据logCPM计算样本相关性和过滤 sample_cor = 'sum,pearson,0'
  if ( ! (sample_cor==FALSE || is.null(sample_cor)) ) { 
    sample_cor <- unlist(str_split(sample_cor,',')) 
    samcor <- SampleCorrelation(counts=normalization$logCPM, sample_group=sample_group, 
                                by='Group', cutoff=as.double(sample_cor[3]), cal=sample_cor[1], method=sample_cor[2], 
                                save=norm_prefix)
    sample_group <- samcor$group
    # correlation_matrix <- samcor$correlation
    # write.table(samcor$count, paste(correlation_prefix,'matrix.tsv',sep='.'), sep="\t", quote=F, col.names=NA, fileEncoding='utf-8')
  }
  
  # --- Condition testing --- #
  
  ClinicalTesting( normalization$DGEListFiltered, sample_group, outdir=outdir, 
                   numeric_cols=numeric_cols, category_cols=category_cols)
  
  # --- Voom correction --- #
  
  norm_exp <- Correction( normalization$DGEListFiltered, sample_group, outdir=outdir,
                          correction_factors=correction_factors)
  
  # --- Differential expression genes --- #
  
  degene_prefix <- OutputPrefix( rootdir=outdir, outdir='degenes', outprefix='degenes' )
  sample_group$Group <- gsub("-", "", sample_group$Group)
  degene_list <- DEGbyCofactor( norm_exp=norm_exp, sample_group=sample_group, cofactors=de_factors, out=degene_prefix )
  
  
  # --- Correlation of DE genes --- #
  group_degene <- degene_list$DEG
  degene_matrix <- norm_exp[rownames(norm_exp) %in% rownames(group_degene),]$E
  feature_correlation <- FeatureCorrelation(degene_matrix, out=paste(degene_prefix,gsub(" ","",degene_list$formula),"sig",sep='.'))

  # 差异基因的聚类效果
  MDSforMatrix(degene_matrix, sample_group[colnames(degene_matrix),], conditions=category_cols, save=paste(degene_prefix,gsub(" ","",degene_list$formula),"sig",sep='.'))
  PCAforMatrix(degene_matrix, sample_group[colnames(degene_matrix),], conditions=category_cols, save=paste(degene_prefix,gsub(" ","",degene_list$formula),"sig",sep='.') )
  tSNEforMatrix(degene_matrix, sample_group[colnames(degene_matrix),], conditions=category_cols, save=paste(degene_prefix,gsub(" ","",degene_list$formula),"sig",sep='.') )
  
  return( list('preprocess'=preprocess, 'normalization'=normalization, 'correction'=norm_exp, 
               'degene'=degene_list, 'degene_cor'=feature_correlation, 
               'group'=sample_group) )
  
}



# ------------------------------------------------------

# START 

# all_cofactors <- c('Group', "Age", "Conservation_days", "RIN_value", "cDNA_conc", "Lib_conc", "Lib_amount", "Data_amount", "Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Sequencing_batch")

outdir <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\limma_RNA_healthy"
count_file_proteincoding <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\input_matrix\\merged.htseq.protein.count.merge_two_batch_used"
count_file_raw <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\input_matrix\\merged.htseq.raw.count.merge_two_batch_used"
group_file <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\input_matrix\\conditions.with_cofactor.all_sample_two_batch.addHealthy.tsv"

Sarcoma_nonSarcoma_oldfilter_proteincoding <- Main ( count_file_proteincoding, 
         group_file, 
         outdir=paste(outdir,'Sarcoma_nonSarcoma_oldfilter_proteincoding',sep="\\"), 
         gene_filter = FALSE, 
         sample_filter = '5,5000', 
         sample_cor = 'sum,pearson,0',
         correction_factors=c('Group', "Age", "RIN_value"), 
         de_factors=c('Group', "Age", "RIN_value", "Batch"),
         numeric_cols=c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'),
         category_cols=c("Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Clinical_group", "Group") )

Sarcoma_nonSarcoma_oldfilter_raw <- Main ( count_file_raw, 
         group_file, 
         outdir=paste(outdir,'Sarcoma_nonSarcoma_oldfilter_raw',sep="\\"), 
         gene_filter = FALSE, 
         sample_filter = '5,5000', 
         sample_cor = 'sum,pearson,0',
         correction_factors=c('Group', "Age", "RIN_value"), 
         de_factors=c('Group', "Age", "RIN_value", "Batch"),
         numeric_cols=c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'),
         category_cols=c("Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Clinical_group", "Group") )

Sarcoma_nonSarcoma_newfilter_raw <- Main ( count_file_raw, 
         group_file, 
         outdir=paste(outdir,'Sarcoma_nonSarcoma_newfilter_raw',sep="\\"), 
         gene_filter = FALSE, 
         sample_filter = '5,5000', 
         sample_cor = 'sum,pearson,0',
         correction_factors=c('Group', 'Conservation_days', 'cDNA_conc', 'Lib_conc', 'Data_amount'), 
         de_factors=c('Group', "Age", "RIN_value", 'Conservation_days', 'cDNA_conc', 'Lib_conc', 'Data_amount', "Batch", "Extraction_batch"),
         numeric_cols=c('Age', 'Conservation_days', 'RIN_value', 'cDNA_conc', 'Lib_conc', 'Lib_amount', 'Data_amount'),
         category_cols=c("Batch", "Hemolysis", "Extraction_batch", "Reverse_batch", "Clinical_group", "Group") )

save.image("limma_RNA_healthy.RData")




# # 条件测试1：查看batch和group两个影响因子交互时差异基因的情况
# new_formula <- "~ 0 + Group*Batch + Age + RIN"
# new_design_model <- model.matrix(as.formula(new_formula), sample_group)
# colnames(new_design_model) <- gsub("Group", "", colnames(new_design_model))
# # > head(new_design_model)
# #           Leiomyoma Sarcoma Batchsecond Age RIN Sarcoma:Batchsecond
# # 990060003         1       0           0  32 1.5                   0
# # 990060015         1       0           0  59 7.3                   0
# # 990060022         1       0           0  76 7.2                   0
# # 990060023         1       0           1  50 7.4                   0
# 
# # Sarcoma:Batchsecond表示Sarcoma和Leiomyoma在second batch和first batch的差异
# # new_contr.matrix <- makeContrasts(Sarcoma - Leiomyoma - Sarcoma:Batchsecond, levels = new_design_model)
# group_batch_degene <- PerformDEanalysis( EList=norm_exp, design_model=new_design_model, coefficients=6, out=paste("count",new_formula,sep="."), pvalue=0.05, logFC=1, plot=F )
# 
# 
# # 条件测试2：只看年龄的模型,未产生差异基因
# Age <- sample_group[colnames(norm_exp),'Age']
# age_model <- model.matrix(~Age)
# rownames(age_model) <- colnames(norm_exp)
# age_degene <- PerformDEanalysis( EList=norm_exp, design_model=age_model, coefficients=2, out=paste("count","Age",sep="."), pvalue=0.05, logFC=1, plot=F )
# 
# 
# # 条件测试3：只看RIN值,未产生差异基因
# RIN <- sample_group[colnames(norm_exp),'RIN']
# rin_model <- model.matrix(~RIN)
# rownames(rin_model) <- colnames(norm_exp)
# rin_degene <- PerformDEanalysis( EList=norm_exp, design_model=rin_model, coefficients=2, out=paste("count","RIN",sep="."), plot=F )
# 
# 
# # 条件测试4：只看batch显著的差异基因
# batch_formula <- "~ 0 + Batch + Age + RIN + Group"
# batch_design_model <- model.matrix(as.formula(batch_formula), sample_group)
# colnames(batch_design_model) <- gsub("Batch", "", colnames(batch_design_model))
# batch_contr.matrix <- makeContrasts(second - first, levels = batch_design_model)
# batch_degene <- PerformDEanalysis( EList=norm_exp, design_model=batch_design_model, contr.matrix=batch_contr.matrix, out=paste("count",batch_formula,sep="."), pvalue=0.05, logFC=1, plot=F )


# 查看样本表达量中位值和年龄、RIN的相关性
# matrix <- norm_exp$E
# 
# library(ggplot2)
# group <- sample_group[colnames(matrix),]
# 
# for ( grp in unique(group$Group) ) {
#   subgroup <- group[group$Group==grp, ]
#   for ( type in c('Age', 'RIN') ) {
#     x <- apply(matrix[,rownames(subgroup)], 2, mean)
#     y <- subgroup[,type]
#     df <- data.frame(x=x, y=y)
#     print(ggplot(df, aes(x = x, y = y)) +
#             geom_point() +
#             stat_smooth(method = "lm", col = "red") +
#             labs(x="Normalized expression", y=type, title = grp)
#     )
#     readline(prompt="Press [enter] to continue")
#   }
# }


# # 条件测试5：对不同分组的数据，如果有NA值则赋值组内平均值
# sample_group <- FillNAvalues( sample_group, by='Group', fill='mean')
# dge_counts <- DGEList(counts, samples=sample_group, genes=gene_length$Length)
# dge_counts$samples <- dge_counts$samples[,-1]
# dge_counts_filtered <- calcNormFactors(dge_counts_filtered, method = "TMM")
# # voom校正数据
# design <- model.matrix(~ 0 + Group + Age + RIN, sample_group)
# colnames(design) <- gsub("Group", "", colnames(design))
# par(mfrow=c(1,1))
# norm_exp <- voom(dge_counts_filtered$counts[,rownames(design)], design, plot = T)
# # write.table(norm_exp$E, "count.limma_normalized.filledNA.tsv", col.names = NA, row.names=T, sep="\t", quote=F)
# # write.table(sample_group[colnames(norm_exp),], "count.limma_normalized.sample_group.filledNA.tsv", col.names = NA, row.names=T, sep="\t", quote=F)
# # 差异基因检测
# cofactors <- c('Group', 'Age', 'RIN', 'Batch')
# design_factor <- colnames(sample_group)[colnames(sample_group) %in% cofactors]
# design_factor <- design_factor[match(cofactors, design_factor)]
# formula <- paste0("~0+", paste0(design_factor, collapse="+"))
# design_model <- model.matrix(as.formula(formula), sample_group)
# colnames(design_model) <- gsub("Group", "", colnames(design_model))
# # 构建比对矩阵
# contr.matrix <- makeContrasts(Sarcoma - Leiomyoma, levels = design_model)
# group_degene <- PerformDEanalysis( EList=norm_exp, design_model=design_model, contr.matrix=contr.matrix, out=paste("count.filledNA",formula,sep="."), pvalue=0.05, logFC=1, plot=F )
# 


# RUVg找到的差异基因和limma-voom找到的差异基因相关性
# zhao <- read.delim("degene_diff\\degene_list.zhaoy.txt", header=F, stringsAsFactors=F, check.names=F)$V1
# zhao <- zhao[zhao %in% rownames(norm_exp)]
# BoxPlotForGene ( norm_exp, sample_group[colnames(norm_exp),], zhao )
# 
# degene <- rownames(group_degene[1:20,])
# degene <- degene[degene %in% rownames(norm_exp)]
# BoxPlotForGene ( norm_exp, sample_group[colnames(norm_exp),], degene )



# 查看重复样本的PCA
# library(FactoMineR)
# library(factoextra)
# repeat_samples <- c("990060033","990060097","990060154","990410002","990410060","990410092","990410096","990410116")
# repeat_samples <- c(paste(repeat_samples, "1", sep="-"), paste(repeat_samples, "2", sep="-"))
# repeat_samples <- repeat_samples[repeat_samples %in% colnames(dge_counts_filtered)]
# subgroup <- sample_group[repeat_samples,]
# subgroup[,c('SampleName', 'SampleReplicate')] <- str_split_fixed(rownames(subgroup), '-', 2)
# new_pca <- PCA(t(dge_counts_filtered$counts[,repeat_samples]), graph = FALSE)
# fviz_pca_ind(new_pca, geom.ind = "text",
#                      col.ind = subgroup[,'SampleName'], addEllipses=T, legend.title='SampleName', title="")



# 分类型数据 Extraction_batch
# library(varhandle)
# extraction_binary <- to.dummy(sample_group[,'Extraction_batch'], "Extraction")
# rownames(extraction_binary) <- rownames(sample_group)
# extraction_test <- data.frame()
# for ( run in colnames(extraction_binary) ) {
#   extraction_binary[,run][extraction_binary[,run]==0] <- "other"
#   extraction_binary[,run][extraction_binary[,run]==1] <- run
#   test <- fisher.test(table(sample_group[, "Group"], extraction_binary[,run]))
#   extraction_test[run,"Pvalue"] <- test$p.value
#   print(table(sample_group[, "Group"], extraction_binary[,run]))
# }
# write.table(extraction_test, paste(test_prefix, "group_class_clinical.fishertest.","Extraction_batch","tsv",sep="."), col.names=NA, sep="\t", quote=F)



