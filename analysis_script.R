install.packages("DatABEL",repos = "http://R-Forge.R-project.org")
BiocManager::install("GenABEL")
rm(list = ls() )
library(rlang)
library(digest)
library(DESeq2)
library(hgu133plus2.db)
library(Biobase)
library(limma)
require(lumi)
library(sva)
library(rafalib)
library(cluster)
library(RColorBrewer)
library(geneplotter)
library(ggfortify)
library(data.table)
library(affycoretools)
library(gplots)
library(gridExtra)
library(ggplotify)
library(tidyr )
library(dplyr)
library(dict)
library(cluster)
library(DatABEL)
library(MetaIntegrator)
expr_m <- clean_data_fin_t

exprs_analysis <- function(expr_m){

  n_cols <- ncol(expr_m) - 2    # Final two columns containing phenotype data
  t_cols <- ncol(expr_m)        # rest of the dataframe to be converted to the expression matrix
  
  #replace sample name to numbers
  dataset_exprs <- expr_m[ ,1:n_cols ]
  dataset_exprs <- apply(dataset_exprs,1,as.numeric)  # Expression matrix should be numeric
  rownames(dataset_exprs) <- colnames(expr_m[,1:n_cols])
  pheno         <- as.data.frame(cbind(replicate(nrow(expr_m),1),expr_m["case"]))  
  pheno$case    <- ifelse(pheno$case == "normal",0,1)
 
  #pheno           <- pheno[,1:2]
  colnames(pheno) <- c("control", "case")
  pheno           <- pheno[order(- pheno$case), ]
  #rownames(dataset_exprs)<-c(1:dim(expr_m)[1])
  # setup experimental desgin
  #design_dataframe  <- as.data.frame(cbind(replicate(nrow(pheno))))
  design            <- model.matrix(~pheno[,2])
  # change colnames to Case and Control
  colnames(design) <- c("intercept", "case")
  # transpose dataset, convert to numeric 
  transposed_dataset_exprs <- dataset_exprs
  #run diff expression
  dataset_exprs_fit             <- lmFit(transposed_dataset_exprs, design)                         # fit linear model for each gene with robust estimation of the variance
  dataset_exprs_contrast_matrix <- makeContrasts(case-intercept, levels=design)                      # create design matrix 
  dataset_exprs_contrast_fit    <- contrasts.fit(dataset_exprs_fit, dataset_exprs_contrast_matrix) # contrasts between the two datasets
  # Extract a table of top ranked genes 
  dataset_exprs_ebayes_fit      <- eBayes(dataset_exprs_contrast_fit, robust = T )

  dataset_exprs_top_genes_wc    <- topTable(dataset_exprs_ebayes_fit,number = nrow(transposed_dataset_exprs), adjust.method = "fdr", confint = TRUE)
  #dataset_exprs_top_genes_int   <- topTreat(dataset_exprs_treat_fit, coef = 2,number = nrow(transposed_dataset_exprs), adjust.method = "fdr", confint = TRUE)
  
   
  return(dataset_exprs_top_genes_wc)
}

expr_m <- clean_data_fin_t
sex_exprs_analysis <- function(expr_m){
  
  n_cols <- ncol(expr_m) - 2    # Final two columns containing phenotype data
  t_cols <- ncol(expr_m)        # rest of the dataframe to be converted to the expression matrix
  
  #replace sample name to numbers
  dataset_exprs <- expr_m[ , 1:n_cols ]
  dataset_exprs <- apply(dataset_exprs,1,as.numeric)  # Expression matrix should be numeric
  rownames(dataset_exprs) <- colnames(expr_m[,1:n_cols])
  
  pheno            <- as.data.frame(cbind(replicate(nrow(expr_m),1),expr_m["case"], expr_m["sex"]))  
  #pheno$case    <- ifelse(pheno$case == "normal",0,1)
  
  #pheno           <- pheno[,1:2]
  colnames(pheno) <- c("control", "case","sex")

  #rownames(dataset_exprs)<-c(1:dim(expr_m)[1])
  # setup experimental desgin
  #design_dataframe  <- as.data.frame(cbind(replicate(nrow(pheno))))
  pheno[,2]        <- factor(pheno[,2], levels = c("normal", "Alzheimer's disease"))
  design            <- model.matrix(~pheno[,2]*pheno[,3] )
  # change colnames to Case and Control
  colnames(design) <- c("intercept", "case", "male","case_t_male")
  # transpose dataset, convert to numeric 
  transposed_dataset_exprs <- dataset_exprs
  #run diff expression
  dataset_exprs_fit             <- lmFit(transposed_dataset_exprs, design)                         # fit linear model for each gene with robust estimation of the variance
# contrasts between the two datasets
  # Extract a table of top ranked genes 
  dataset_exprs_ebayes_fit      <- eBayes(dataset_exprs_fit, robust = T )
  #contrast_fit_ebayes_int       <- contrasts.fit(dataset_exprs_ebayes_fit, c(0,0,-1,1))
  #final_fit_interaction         <- eBayes(contrast_fit_ebayes_int)
  dataset_exprs_top_genes_wc    <- topTable(dataset_exprs_ebayes_fit,coef =4, adjust.method = "fdr", confint = TRUE,number = nrow(transposed_dataset_exprs))
  
  
  
  return(dataset_exprs_top_genes_wc)
}

dataset <- expr_ae_ad
extract_good_probe_list <- function(dataset) {
  # dataset - expression dataset as dataframe
  # probe_percentile_threshold - percentile at which to use as cut-off for detected probes 
  # number of samples in which probe must be expressed in - fixed at 0.8 - i.e 80% of samples
  # calculate quantile threshold for each sample
   # samples as columns and genes as rows
  dataset <- t(dataset)
  sample_quantiles <- apply(dataset, 2, quantile, probs=0.90 ) # Find the 90th percentile of for each sample
  sample_quantiles <- replicate(nrow(dataset), sample_quantiles)
  sample_quantiles <- t(sample_quantiles)
  # convert to dataframe
  dataset_check     <- dataset- sample_quantiles
  dataset_check     <- as.data.frame(dataset_check)
  dataset_replace   <- dataset_check %>% mutate_each(funs(replace(.,.<0,NA)))
  
  rownames(dataset_replace) <- rownames(dataset)
  count_na_row              <- apply(dataset_replace,1,function(x){sum(is.na(x))})
  dataset_check_2           <- as.data.frame(cbind(dataset, count_na_row))
  dataset_check_2$pct       <- dataset_check_2$count_na_row > 0.8 * ncol(dataset)
  dataset_check_2           <- subset(dataset_check_2, dataset_check_2$pct == FALSE)

  # subset good probes
  good_probes <- rownames( dataset_check_2)

  return(good_probes)
}


sva_cleaning_mat <- function(edat){
  # edat is the expression matrix with sex and case as the last two columns 
  mod_data   <- model.matrix( ~1 + case + sex, data = edat)
  mod0       <- model.matrix(~1+ sex, data = edat)
  n_cols     <- ncol(edat) - 2
  p_cols     <- ncol(edat) - 1
  expr_ae_ad <- apply(edat[,1:n_cols],2,as.numeric)
  #expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
  #expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
  
  n.survar   <- num.sv( t(expr_ae_ad),mod_data, method = "leek") 
  if(n.survar > 0){
      svobj  <- sva(t(expr_ae_ad), mod_data, mod0, n.sv=n.survar, method="two-step")
  # adjust for sva
      X     <- cbind(mod_data, svobj$sv)
      Hat   <- solve(t(X) %*% X) %*% t(X)
      beta  <- (Hat %*% expr_ae_ad)
      P     <-  ncol(mod_data)
  
      clean_data <- t(expr_ae_ad) - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),])
      clean_data <- t(clean_data)
      clean_data <- as.data.frame(cbind(clean_data,as.character(edat[,p_cols]),as.character(edat[,ncol(edat)])))
      colnames(clean_data) <- colnames(edat)
      rownames(clean_data) <- rownames(edat)
  }
  else{
      
      clean_data <- edat
    
  }
  return(clean_data)
} 



#source ("M:/Consulting_Core/Pathway_methylation/DATASETS/ROSMAP/plotPCA.R")
# Robust Spline Normalization, Background Correction and Log2 transformation of the expression matrix, use IF AND ONLY IF
# you are unable to extract expression data from "good_run_script.R"

 n_m_two    <- ncol(eset_ae_ad_hip) - 2
#expr_ae_ad  <- eset_ae_ad_p[,1:n_m_two]

expr_ae_ad <- as.matrix(eset_ae_ad_hip[,1:n_m_two])
colnames(expr_ae_ad) <- colnames(eset_ae_ad_hip[,1:n_m_two])
rownames(expr_ae_ad) <- rownames(eset_ae_ad_hip[,1:n_m_two])
expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)
#rownames(expr_ae_ad) <- rownames(eset_ae_ad_t)ls
#colnames(expr_ae_ad) <- colnames(eset_ae_ad_t)

probes_names  <- extract_good_probe_list(t(expr_ae_ad))
rownames(eset_ae_ad) <- colnames(expr_ae_ad)
colnames(eset_ae_ad) <- rownames(expr_ae_ad)
rowname_vec2  <- rownames(eset_ae_ad)

expr_ae_ad    <- expr_ae_ad[,colnames(expr_ae_ad) %in% probes_names]


eset_ae_ad_pb <- as.data.frame(cbind(expr_ae_ad,hip_pData$FactorValue.disease.,hip_pData$Characteristics..sex.))
colnames(eset_ae_ad_pb) <- c(probes_names, "case", "sex")
#colnames(eset_ae_ad_pb)  <-  AE_AD_pData$IID


rownames(eset_ae_ad_pb) <- hip_pData$Source.Name

clean_mat_ad    <- sva_cleaning_mat(eset_ae_ad_pb)
sva_check       <- sva_cleaning_mat(clean_mat_ad)
#expr_ae_ad <- as.matrix(eset_ae_ad_pb[,1:2063])
#expr_ae_ad <- apply(expr_ae_ad,2,as.numeric)


#clean_data_m           <- as.data.frame(cbind(t(clean_mat_ad)))
clean_data             <- as.data.frame(clean_mat_ad)
clean_data_m           <- as.data.frame(cbind(clean_mat_ad, hip_pData$FactorValue.disease., hip_pData$Characteristics..sex.))

colnames(clean_data_m) <- c(colnames(clean_mat_ad), "case","sex")



# matrix with batch effect corrections 



# Need to annotate chip specific names to Entrez IDs


exp_df_wo_s  <- exprs_analysis(clean_data_final)[[1]]
#exp_df_wo_s  <- subset(exp_df_wo_s, exp_df_wo_s$adj.P.Val < 0.05 )
clean_data_final <- as.data.frame(t(as.matrix(clean_data_fin)))
c_names          <- colnames(clean_data_final)
clean_data_final <- as.data.frame(cbind(clean_data_final, AE_AD_pData$FactorValue..disease., AE_AD_pData$Characteristics.sex.))

colnames(clean_data_final)  <- c(c_names,"case", "sex")


clean_data_fin_t$case <- ifelse(clean_data_fin_t$case == "AD","AD","normal")
exp_df_male  <- subset(clean_data_fin_t, clean_data_fin_t$sex == "male")
exp_df_fem   <- subset(clean_data_fin_t, clean_data_fin_t$sex == "female")

results_limma_male        <- exprs_analysis(exp_df_male)
results_limma_female      <- exprs_analysis(exp_df_fem)

results_limma_male_comp  <- exprs_analysis(exp_df_male)
results_limma_male_t     <- exprs_analysis(exp_df_male)[[2]]


exp_df_ma$fc      <- 2^exp_df_ma$logFC
exp_df_f$fc      <- 2^exp_df_f$logFC
exp_df_wo_s$fc   <- 2^exp_df_wo_s$logFC 
expression_matrix <- clean_data_final[,1:2742]
expression_matrix <- apply(expression_matrix,2,as.numeric)
sst <- rowSums(t(expression_matrix)^2)
ssr <- sst - results_limma_male_comp[[2]]$df.residual * results_limma_male_comp[[2]]$sigma^2
file_name
rsq <- ssr/sst

write.csv(results_limma_male, paste0("male_hip_",file_name,".csv"))
write.csv(results_limma_female, paste0("female_hip_",file_name,".csv"))
save.image(paste0(file_name, "_hip.RData"))

write.table(x = as.matrix(clean_data_fin_t),file = paste0(file_name,"_eset_final.txt"), sep = "\t",row.names = TRUE, col.names = TRUE)
