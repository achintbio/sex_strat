# Script changes machine specific probe names to ENTREZ Gene ids. Multiple probes mapped to the same entrez ID were picked based on larger variance across
# samples, need to input the type of chip sequencing technology used the and the initial dataset
BiocManager::install("illuminaHumanWGDASLv4.db")
BiocManager::install("GO.db")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("illuminaHumanv3.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("hgu133a.db")
BiocManager::install("illuminaHumanv4.db")
BiocManager::install("hugene10sttranscriptcluster.db")

library(annotate)
library(AnnotationDbi)
library(illuminaHumanWGDASLv4.db)
library(hgu133plus2.db)
library(hugene10stv1cdf)
library(hugene10sttranscriptcluster.db)
library(illuminaHumanv3.db)
library(GO.db)
library(org.Hs.eg.db)
library(pd.hugene.1.0.st.v1)

largest_var_picker <- function(x){
  dataset_expr_c <- subset(clean_data_t, clean_data_t$entidc == x) 
  dataset_expr_c <- dataset_expr_c[ order(dataset_expr_c$var), ]
  dataset_expr_c$id <- rownames(dataset_expr_c)
  return(dataset_expr_c[1,ncol(dataset_expr_c)]) 
}

clean_data_fin_t <- clean_data_m
annotation_object <- hugene10sttranscriptclusterSYMBOL

chip_entrez  <- hugene10sttranscriptclusterENTREZID
annotation_function <- function(clean_data_fin_t,annotation_object) { 
  
  
  
   
   annotation_object_ent  <- annotation_object                       # Annotation Object converting hgu133plus2 specific probe IDS to Entrez Gene Names
   mapped_probes          <- mappedkeys(annotation_object_ent)                 # find the probe IDS which have a mapping to the Entrez ID   
#mapped_probes_gol     <- as.list(annotation_object_sym[mapped_probes])


  clean_data_t          <- as.data.frame(t(as.matrix(clean_data_fin_t)))       # convert expression matrix to a dataframe with the genes as rows 
  clean_data_t$chk      <- rownames(clean_data_t) %in% mapped_probes        # find which genes are mapped 
  clean_data_t          <- subset(clean_data_t, clean_data_t$chk == TRUE)   # drop genes which aren't mapped

  clean_data_t$entid   <- mapped_probes_l[rownames(clean_data_t)]           # Find the entrez IDS
  clean_data_t$entidc  <- as.character(clean_data_t$entid)                  # Convert to character         # 


  expr_mat_clean      <- as.matrix(clean_data_t[,1:nrow(hip_pData)])                     # Drop the redundant columns
  var_by_sample       <- apply(expr_mat_clean,1, var)                        # Find Variance by row to resolve ties between probes
  clean_data_t$var    <- var_by_sample                                       # Add variances to the expression matrix
  clean_data_t$dup    <- duplicated(clean_data_t$entidc)                     # Find which entrez IDs are duplicated
  dup_ID_df           <- subset(clean_data_t, clean_data_t$dup == TRUE)      # create a data frame with all duplicated ids
  ent_id_occ          <- dup_ID_df$entidc                                    # 
  clean_data_t$dup    <- clean_data_t$entidc %in% ent_id_occ
  clean_data_t        <- clean_data_t[ order(clean_data_t$var), ]


  ill_id_lar_var     <- lapply(ent_id_occ, largest_var_picker)               # Resolve ties by picking the probe with larger variance
  ill_id_lar_var     <- as.character(ill_id_lar_var)                         #    
  clean_data_dup     <- subset(clean_data_t, clean_data_t$dup == TRUE)
  clean_data_dup$id  <- rownames(clean_data_dup)                             # Create a column to check if the datasets have
  clean_data_dup$lv  <- clean_data_dup$id %in% ill_id_lar_var                # Check if the duplicated IDs are the same with the larger vairance  
  clean_data_dup     <- subset(clean_data_dup, clean_data_dup$lv == TRUE)    # Keep the IDS with the larger Variance  
  clean_data_ndup    <- subset(clean_data_t, clean_data_t$dup == FALSE)
  names_row          <- rbind(matrix(clean_data_dup$entidc,ncol =1), matrix(clean_data_ndup$entidc,ncol = 1))
  clean_data_dup     <- clean_data_dup[,1:nrow(hip_pData)]
  clean_data_ndup    <- clean_data_ndup[,1:nrow(hip_pData)]

  clean_data_fin           <- rbind(clean_data_dup, clean_data_ndup)
  rownames(clean_data_fin) <- names_row
#clean_data_fin           <- clean_data_fin[,1:ncol(eset_AE_AD)]

  clean_data_fin_t             <- as.data.frame(cbind(t(as.matrix(clean_data_fin)),as.character(hip_pData$FactorValue.disease.),as.character( hip_pData$Characteristics..sex.)))
  colnames(clean_data_fin_t)   <- c(rownames(clean_data_fin), "case" ,"sex") 
  rownames(clean_data_fin_t)   <- hip_pData$Source.Name

return(clean_data_fin_t)
} 
x <- colnames(clean_data_fin_t)


clean_data_fin_t <- eset_ae_ad_tx
clean_data_fin_t <- annotation_function(clean_data_m, hugene10sttranscriptclusterSYMBOL)
