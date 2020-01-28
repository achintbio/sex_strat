##########################################################################################################################################
####                                                                                                                                  ####
###                                                                                                                                    ###
##                                                                                                                                      ##
#                                                  run_script_loaded_data_file                                                            #
##                                                                                                                                      ##
###                                                                                                                                    ###
####                                                                                                                                  ####
##########################################################################################################################################
# AUTHOR: Achintya VARMA
# QC PIPELINE VERSION: 2.0
# DATE: 11/08/2019
# ARRAY EXPRESS NUMBER: E-GEOD-36980
# DISORDER: Alzheimer's Disease
# MICROARRAY PLATFORM: Affymetrix
# EXPRESSION CHIP: Human Genome U133A
# NUMBER OF SAMPLES: 79
# TISSUE: Hippocampus, Temporal Cortex and Frontal Cortex
#
# NOTES - 
# Read RAW DATA, extract phenotypes and find data related to the study
# e.g given using E-GEOD-3970 dataset
# 

##### SET PARAMETERS #####
rm( list = ls() )
library(ArrayExpress)
require(affycoretools)
library(rlang)
library(Biobase)
library(affy)
library(oligo)
library(devtools)
library(massiR)
library(annotate)
require(metap)
require(lumi)
library(hugene10sttranscriptcluster.db)
library(hgu133plus2.db)
library(pd.hugene.1.0.st.v1)
library(illuminaHumanv4.db)

AE_AD_eset3    <- ae2bioc(mageFiles = raw_data_3) # Read raw data files
exprs_conv     <- oligo::rma(AE_AD_eset3)         # Robust Multichip averaging -> returns background corrected, normalized and log2 expression converted expression sets


eset_ae_ad     <- exprs(exprs_conv)
row_name_v     <- rownames(eset_ae_ad)
eset_ae_ad     <- apply(eset_ae_ad,2,as.numeric)

rownames(eset_ae_ad)  <-  row_name_v
#eset_ae_ad     <- apply(eset_ae_ad,1,log2)
AE_AD_pData    <- pData(phenoData(exprs_conv))           # get phenotype data
AE_AD_pData    <- subset(AE_AD_pData, 
                         AE_AD_pData$FactorValue..organism.part. %in% 
                         c("hippocampus","entorhinal cortex","temporal cortex", "middle temporal gyrus"))            # subset to include only tissues from temporal lobe
 

rownames(AE_AD_pData)   <- gsub(".CEL", "" , rownames(AE_AD_pData))           
colnames(eset_ae_ad)    <- gsub(".CEL", "" , colnames(eset_ae_ad))

eset_ae_ad_t            <- as.data.frame(t(eset_ae_ad))


eset_ae_ad_t$link_expr_pheno <- rownames(eset_ae_ad_t) %in% rownames(AE_AD_pData)
eset_ae_ad_t            <- subset(eset_ae_ad_t, eset_ae_ad_t$link_expr_pheno == TRUE)
eset_ae_ad_t            <- eset_ae_ad_t[ ,1:nrow(eset_ae_ad)]

rm(eset_ae_ad_f)

eset_ae_ad_p            <- as.data.frame(cbind(eset_ae_ad_t, AE_AD_pData$Characteristics..sex. ,AE_AD_pData$Factor.Value..disease.)) # combine gene expression matrix with case status and gender

colnames(eset_ae_ad_p)  <- c(colnames(eset_ae_ad_t), "sex", "disease" )






save.image(paste0(file_name,".","RData"))


