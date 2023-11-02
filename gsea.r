# GSEA两个基因列表富集程度
# ref: https://cloud.tencent.com/developer/article/1838918
#      https://www.jianshu.com/p/cb5afe6aed70
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(enrichplot)

genelist_df_file <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\limma_RNA_healthy\\Sarcoma_nonSarcoma_oldfilter_proteincoding\\degenes\\degenes.~0+Group+Age+RIN_value+Batch.tsv"
# genelist_query_file <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2\\分类错误样本比较\\gene_list.zhao.tsv"
genelist_query_file <- "F:\\项目\\TEP子宫肉瘤\\tpm_matrix\\merge_batch1_batch2_addHealthy\\limma_RNA_healthy\\Sarcoma_nonSarcoma_oldfilter_proteincoding\\degenes\\degenes.~0+Group+Age+RIN_value+Batch.sig.tsv"

genelist_db <- read.table(genelist_df_file, sep="\t",header=T, row.names = 1, check.names = F)
genelist_db.id <- bitr(rownames(genelist_db), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dup_id <- genelist_db.id$SYMBOL[duplicated(genelist_db.id$SYMBOL)]

genelist_db <- genelist_db[!rownames(genelist_db) %in% dup_id,]
genelist_db.id <- genelist_db.id[!rownames(genelist_db) %in% dup_id,]
rownames(genelist_db.id) <- genelist_db.id$SYMBOL
genelist_db.id$logFC <- genelist_db[match(rownames(genelist_db), genelist_db.id$SYMBOL), 'logFC']

genelist_db.genelist <- genelist_db.id$logFC
names(genelist_db.genelist) <- genelist_db.id$ENTREZID
genelist_db.genelist <- sort(genelist_db.genelist, decreasing = T)

# genelist_query <- read.table(genelist_query_file, header=F)$V1
genelist_query <- read.table(genelist_query_file, sep="\t",header=T, row.names = 1, check.names = F)
genelist_query <- bitr(rownames(genelist_query), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
genelist_query.genelist <- seq(-1,-length(genelist_query$ENTREZID))
names(genelist_query.genelist) <- genelist_query$ENTREZID

# 两个基因列表的比较
genelist_query2db <- GSEA(genelist_db.genelist, TERM2GENE=data.frame('term'='query.genelist', 'gene'=names(genelist_query.genelist)), 
                        pvalueCutoff=1, minGSSize = 3 ) #GSEA
gseaplot2(genelist_query2db, geneSetID=1, title=paste0("DB genelist (N=",length(genelist_db.genelist),")"), pvalue_table=T)


# gsea
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
GeneSet <- read.gmt(gmtfile) ##读取gmt文件得到基因集S
genelist_query.gsea <- GSEA(genelist_query.genelist, TERM2GENE=GeneSet, pvalueCutoff=0.5 ) #GSEA
head(genelist_query.gsea)
dotplot(genelist_query.gsea)

# kegg
gse.KEGG <- gseKEGG(genelist_query.genelist, 
                    organism = "hsa", # 人 hsa
                    keyType = "kegg",
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH") 
gseaplot2(gse.KEGG, 1)

# go
genelist_query.go <- gseGO(
  genelist_query.genelist, #geneList
  ont = "BP",  # 可选"BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db, #人 注释基因
  keyType = "ENTREZID",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",#p值校正方法
  minGSSize =10, maxGSSize=1000
)
gseaplot2(genelist_query.go, geneSetID=1:5, title=genelist_query.go$Description[1:5], pvalue_table=T)
ridgeplot(genelist_query.go)

