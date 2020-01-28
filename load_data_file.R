##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                  DOWNLOAD E-GEOD- 36980                                                                 #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################
# AUTHOR: Achintya VARMA
# QC PIPELINE VERSION: 1.0
# DATE: 10/07/2019
# ARRAY EXPRESS NUMBER: E-GEOD-36980
# DISORDER: Alzheimer's Disease
# MICROARRAY PLATFORM: Affymetrix
# EXPRESSION CHIP: ?
# NUMBER OF SAMPLES: 79
# TISSUE: Hippocampus, Temporal Cortex and Frontal Cortex
#
# NOTES - 
# Read RAW DATA, extract phenotypes and find data related to the study
# e.g given using E-GEOD-3970 dataset
# 

##### SET PARAMETERS #####


rm(list = ls())
library(beadarray)
library(annotat)
library(AnnotationDbi)
require(ArrayExpress)
require(affy)
require(icesTAF)
library(biomaRt)
library(arrayQualityMetrics)
library(Rcpp)
library(limma)
library(sva)
library(statmod)
library(readr)
library(readxl)
library(synapser)
library(synapserutils)
synLogin(email = "amvarma123", password = "ManBaapJi&2018")
entity <- synGet("syn3157225",downloadLocation = getwd())

#results <- synTableQuery(sprintf("select * from %s", file_name))
# Need this to make directories in R environment 
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(mart = mart,attributes = c(
 #   "affy_hg_u133_plus_2",
  #  "ensembl_gene_id",
   ##"external_gene_name"),
  #filter = "affy_hg_u133_plus_2",
  #values = rownames(exprs(AE_AD_eset3)),
  #uniqueRows=TRUE)
#mkdir(work_dir)

rm(list = ls())

work_dir          <- "C:/Users/axv851/Documents/Meta_Analysis_Pre_Processing/Analysis/r_images"
#change directory address if needed
setwd(work_dir)
 
file_name         <- "E-GEOD-48350" # Need this ID to download files from ARRAY Express
file_name_raw     <- paste(file_name, "raw.1", sep = ".")

#setwd(paste(work_dir,file_name, sep = "/"))
#mkdir(paste(work_dir,file_name_raw, sep ="/"))
#setwd(paste(work_dir, file_name_raw, sep="/"))
#load(paste0(file_name,".Rdata"))

load(paste0(file_name,".RData"))
sex_inter_effect <- sex_exprs_analysis(clean_data_fin_t)
setwd("C:/Users/axv851/Documents/Meta_Analysis_Pre_Processing/Analysis/interaction_analysis")
write.csv(sex_inter_effect,paste0(file_name, "_interaction.csv"))

raw_data_3                                                  <- getAE(file_name,type = "full")
affy_batch                                                  <- ReadAffy()
#
ad_pdata_print <- as.data.frame(cbind(AE_AD_pData$Source.Name, AE_AD_pData$FactorValue..disease., AE_AD_pData$Characteristics.sex.))
colnames(ad_pdata_print) <- c("ID", "case", "sex")
#AE_AD_pData <- read_excel(paste0(getwd(),"/","mayoGWAS.xls"))
#eset_AE_AD  <- read_csv(paste0(getwd(),"/","mayoEGWAS.csv"))
eset_ae_ad_tx <- read.csv("MayoEGWAS_arrayExpression_TCX.csv")
AE_AD_pData   <- read_xls("mayoGWAS.xls")
link             <- AE_AD_pData$IID %in% eset_ae_ad_tx$IID
AE_AD_pData$link <- link 
AE_AD_pData      <- subset(AE_AD_pData,AE_AD_pData$link == TRUE)
eset_ae_ad_tcx_2       <- subset(eset_ae_ad_tx, eset_ae_ad_tx$IID %in% AE_AD_pData$IID)
eset_ae_ad_tx    <- eset_ae_ad_tcx_2
rm(eset_ae_ad_tcx)
rm(eset_ae_ad_tcx_2)
View(eset_ae_ad_tx)
rownames(eset_ae_ad_tx) <- eset_ae_ad_tx$IID
eset_ae_ad_tx       <- eset_ae_ad_tx[,3:ncol(eset_ae_ad_tx)]
rownames(eset_ae_ad_tx) <- eset_ae_ad_tx$IID
eset_ae_ad_tx       <- as.data.frame(t(as.matrix(eset_ae_ad_tx)))
eset_AE_AD       <- eset_AE_AD[2:nrow(eset_AE_AD), ]

rm(mayoEGWAS)
rm(mayoGWAS)
pheno_29378 <- as.data.frame(cbind(AE_AD_pData$Source.Name, AE_AD_pData$Factor.Value.disease.,AE_AD_pData$Characteristics.sex.))

colnames(pheno_29378) <- c("source_name", "case", "sex")

write.csv(ad_pdata_print, "pheno_28146.csv")
#unique(AE_AD_pData$Comment..Sample_characteristics.)
#AE_AD_eset3                                                 <- ae2bioc(mageFiles = raw_data_3)
#AE_AD_eset2.1                                              <- AE_AD_eset2[[1]]
#unique(AE_AD_pData$FactorValue..POST.MORTEM.DELAY.)