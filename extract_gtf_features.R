#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author  : zhangsiwen
# @ChangeLog
#     20230316, first version


suppressMessages(library(optparse))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(Repitools))
suppressMessages(library(dplyr))
# suppressMessages(library(org.Hs.eg.db))


option_list = list(
  make_option("--gtf", type="character", default=NULL, help="gtf file for feature extraction", metavar="character"),
  make_option("--out", type="character", default=NA, help="out bed file, default {gtf}.bed", metavar="character"),
  make_option("--feature", type="character", default='exon', help="feature, choose from exon | intron | transcript", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


Main <- function ( gtf_file, out=NA, feature='exon' ) {
  if ( is.na(out) ) {
    out <- gsub('gtf', paste(feature,'bed',sep='.'), gtf_file)
  } else {
    outdir <- dirname(out)
    if ( ! dir.exists(outdir) ) {
      dir.create(file.path(outdir), recursive = TRUE)
    }
  }
  
  txdb <- makeTxDbFromGFF(gtf_file, organism='Homo sapiens', circ_seqs='MT')
  
  if ( feature == 'exon' ) {
    extraction <- tidyExons(txdb, drop.geneless=FALSE)
  } else if  ( feature == 'intron' ) {
    extraction <- tidyIntrons(txdb, drop.geneless=FALSE)
  } else if  ( feature == 'transcript' ) {
    extraction <- tidyTranscripts(txdb, drop.geneless=FALSE)
  }
  
  df <- annoGR2DF(extraction)
  
  # add intron number
  if ( feature == 'intron' ) {
    for ( i in 1:nrow(df) ) {
      tx <- df[i,'tx_name']
      df[i,'intron_num'] <- sum(df[1:i,'tx_name'] == tx)
    }
  }
  
  # add gene symbol
  # annot <- select(org.Hs.eg.db,
  #                 keys = keys(org.Hs.eg.db),
  #                 columns = c('ENTREZID','SYMBOL','ENSEMBL','ENSEMBLTRANS'),
  #                 keytype = 'ENTREZID')
  # 
  # geneid <- gsub("[.][0-9*]_[0-9*]$", "", df$gene_id)
  # df <- data.frame(df, 'gene_name'=annot[match(geneid,annot$ENSEMBL), 'SYMBOL'])
  
  # output
  colnames(df)[1] <- paste0('#',colnames(df)[1])
  df$start <- df$start - 1 
  df <- df %>% mutate_if(is.character, trimws) # remove all whitespaces
  df <- format(df, scientific = FALSE)
  output.file <- file(out, "wb")
  write.table(df, file=output.file, quote=F, sep='\t', row.names=F, col.names=T, eol='\n')
  close(output.file)

}


Main ( opt$gtf, out=opt$out, feature=opt$feature )
