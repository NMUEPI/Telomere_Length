#TF analysis

*Based on the single-cell RNA-seq data, R package named Scenic (https://github.com/aertslab/SCENIC) was used to perform transcript factor analysis.*

```R
library(readxl)
gene_symbol<-read_xlsx("GSE136714_raw.counts.for.geo.xlsx")[,1]

library(biomaRt)
library(data.table)
mart <- useMart("ensembl","mmusculus_gene_ensembl")##小鼠选择mmusculus_gene_ensembl
gene_id<-getBM(attributes=c("external_gene_name","ensembl_gene_id"),filters = "ensembl_gene_id",values = gene_symbol, mart = mart)#将输入的filters设置未external_gene_name(也就是gene_symbol),将输出的attributes设置为external_gene_name和emsembl_gene_id
geoData <- fread("8cell_32cell_counts.csv",h=T)

geoData <- read_xlsx("GSE136714_raw.counts.for.geo.xlsx")
geoData1<-merge(geoData,gene_id,by.x="gene",by.y="ensembl_gene_id")

geoData <- read.csv("GSE136714_raw.counts.for.geo_c.csv",stringsAsFactors = F)
geoData<-geoData[c(grep("8cell|32cell|external_gene_name",colnames(geoData)))]

geneNames <- unname(geoData[,138])
exprMatrix <- as.matrix(geoData[,-c(138)])
rm(geoData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames
exprMatrix <- exprMatrix[!duplicated(rownames(exprMatrix)),]
exprMatrix[1:5,1:4]

org<-"mgi"
dbDir<-"/public/user/zj2020/scRNAseq/ref"
ntitle<-"mouse_emryo"
library(SCENIC)
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=ntitle, nCores=10) 

genesKept <- geneFiltering(exprMatrix, scenicOptions=scenicOptions,minCountsPerGene=3*.01*ncol(exprMatrix),minSamples=ncol(exprMatrix)*.01)
exprMat_filtered <- exprMatrix[genesKept,]


interestingGenes <- c("Tert")
interestingGenes[which(!interestingGenes %in% genesKept)]

corrMat <- cor(t(exprMat_filtered), method="spearman")

saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))



logexprMat_filtered <- as.matrix(log2(exprMat_filtered+1))

runGenie3(logexprMat_filtered, scenicOptions)

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

library(BiocParallel)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, as.matrix(log2(exprMat_filtered+1)))
runSCENIC_4_aucell_binarize(scenicOptions)

```