# The following code contains functions for dendogram plotting and pca plots
# Dendogram plots are generated from aggolmerative hiearchical clustering
# The division of datasets based on color in two completely different branches of the cluster signals batch effects
# Outliers in dendograms are classified in their cluster higher up the y axis(i.e - High value of 1 -ISA, where ISA stands for Inter-Sample Adjacency )

# The second function plots out PCA plot- Descriptive plots which aid in assessment of batch effects and outliers. Outliers were defined as points falling out of
# the third standard deviation of the first two principal components

library(ggplot2)
library(ggdendro)
library(WGCNA)
library(DupChecker)

dendogram_plot <- function(dataset) {
  # dataset should be formatted such that the final two columns contain phenotype info- i.e case status and sex 
  # the columns before contain the expression matrix with features as columns and samples as rows
  n_cols  <- ncol(dataset) - 2 # final two columns contain phenotype 

  datExprs  <- t(as.matrix(dataset[,1:n_cols])) # Extract the expression matrix from the dataframe
  datExprs  <- apply(datExprs,2,as.numeric)     # convert the matrix into numeric type matrix
  diagnosis <- dataset[ncol(dataset) - 1]       # extract case status
  sex       <- dataset[ncol(dataset)]           # extract sex status 
  colVec    <-  labels2colors(as.factor(diagnosis)) # Color all the diagnosis into turqoise and blue
  local(
    {colLab <<- function(n,treeorder) {
      if(is.leaf(n)) {
        a <- attributes(n)
        i <<- i+1
        attr(n, "nodePar") <- c(a$nodePar, list(lab.col = colVec[treeorder][i], lab.font = i%%3))
      }
      n
    }
    i <- 0
    })
  # colLab adds colors to the leaves of the dendogram object depending on whether the samples lie in Cases or Controls
 
  cat(" setting up data for qc plots","\r","\n") 
  ## expression matrix and IAC
  cat(" expression matrix and IAC","\r","\n")
  IAC   <- cor(datExprs)   # Find corellations two at a time
  IAC_d <- 1-IAC           # 1 - Corellations
  
  sample_names <- colnames(datExprs)
  IAC        <- cor(datExprs, method="p",use="p")
  diag(IAC)  <- 0
  A.IAC      <- ((1+IAC)/2)^2
  
  cluster1      <- hclust(as.dist(1-A.IAC),method="average") # Use average method to build the dendogram of samples based on expression
  cluster1order <- cluster1$order                            # order of clustering starting from individual samples  
  cluster2      <- as.dendrogram(cluster1,hang=0.1)          # Create a dendogram based on the expression values
  cluster3      <- dendrapply(cluster2,colLab,cluster1order) # color the leaves of the dendogram
  meanIAC       <- apply(IAC,2,mean)                         # find the mean of intersample corellations used for standardization
  sdCorr        <- sd(meanIAC)                               # find the standard deviation : used for standardization
  numbersd      <- (meanIAC - mean(meanIAC))/sdCorr
  numbersd      <- as.data.frame(cbind(rownames(numbersd), numbersd))
  numbersd      <- as.data.frame(cbind(rownames(numbersd), numbersd))
  ind           <- 1:nrow(numbersd)
  numbersd      <- as.data.frame(cbind(numbersd,ind))
  colnames(numbersd) <- c("source_name", "std_corr","ind")
  numbersd      <- numbersd[order(numbersd$source_name), ]
   

  plot(cluster3,nodePar=list(lab.cex=1,pch=NA),
      main=paste("Mean ISA = ",signif(mean(A.IAC[upper.tri(A.IAC)]),3),sep=""),
       xlab="",ylab="1 - ISA",sub="",cex.main=1.8,cex.lab=1.4)
  mtext(paste("distance: 1 - ISA ",sep=""),cex=0.8,line=0.2)
  
 
  # Standardized connectivity plot : plot of standardized connectivity of individual samples vs the index number 
  std_plot <- ggplot(numbersd,aes(x = ind, y = std_corr)) + geom_point() + geom_text(aes(label = numbersd$source_name)) + geom_hline( yintercept = -2) 
  return(list(std_plot, cluster3))
  
}

pca_data_out <- function(expr_m, exp_var){
  
  # function takes a dataframe with expression matrix and phenotype information combined as the first argument and "coloring" phenotype
  # as the second to return a list containing a PCA plot, pca object, list of the four values for outlier detection along the first two components
  # and a dataframe consisting of the first two components and case status. 
  n_cols            <- ncol(expr_m) - 2 # Final two columns contain the expression matrix
  
  gene_exp_mat  <- t(as.matrix(expr_m[,1:n_cols]))  # gene expression matrix
  gene_exp_mat  <- apply(gene_exp_mat,2,as.numeric) # convert the values to numeric 
  gene_exp_mat  <- t(gene_exp_mat)                  # samples as columns and features as rows  
  pca           <- prcomp( gene_exp_mat, center = TRUE, scale. = TRUE ) # build pca object
  pca_df        <- as.data.frame(cbind(pca$x[,1],pca$x[,2],expr_m[exp_var])) # create a dataframe with the first two principal components and the phenotype variable
  colnames(pca_df) <- c("PC1", "PC2","status") 
  pca_df$status <- as.factor(pca_df$status)
  pca_df$pat    <- rownames(pca_df)
   # create a ggplot object with the phenotype column as the color 
  hor1          <- sd(pca$x[,1])
  hor2          <- sd(pca$x[,2])
  ver1          <- mean(pca$x[,1])
  ver2          <- mean(pca$x[,2])
  exp_pca       <- ggplot(pca_df, aes(x = PC1,y = PC2)) + geom_point(aes(color = status)) 
  exp_pca       <- exp_pca + geom_vline( xintercept = ver1 + 3*hor1 ) + geom_vline( xintercept = ver1 - 3*hor1)    # add vertical lines to cover +- 3 standard deviations
  exp_pca       <- exp_pca + geom_hline( yintercept = ver2 + 3*hor2)  +  geom_hline( yintercept = ver2 - 3 *hor2) # add horizontal lines to cover +-3 standard deviations
  #exp_pca       <- exp_pca + geom_label(data = filter(pca_df,pca_df$PC1 > ver1 + 3*hor1), aes(label = pca_df$pat))
  #ht_map        <- coolmap(t(gene_exp_mat),by = "expression level", col = 'redgreen', show.dendrogram = "column",main = paste0("heatmap of ",exp_var))
  #ht_map        <- heatmap.2(gene_exp_mat, Rowv = TRUE, distfun = dist, hclustfun = hclust, dendogram = c("column"), scale = "column", col = my_pallette)
  #   ,  main = paste0("heatmap of",exp_var),  col = "heat.colors")
  
  
  
  return(list(exp_pca, pca,list(hor1,hor2,ver1,ver2),pca_df))
  
}

clean_data_pre_ann  <- as.data.frame(cbind(clean_mat_ad,AE_AD_pData$Factor.Value.disease., AE_AD_pData$Characteristics.sex.))
colnames(clean_data_pre_ann)  <- c(probes_names, "case", "sex")
pca_out_c  <- pca_data_out(clean_data_m, "case")
pca_out_s  <- pca_data_out(clean_data_m,  "sex")
pca_out_s
 

clean_data_m <- clean_data_m[!(rownames(clean_data_m) %in% "GSM907837 1"), ]
clean_data_m <- clean_data_m[!(rownames(clean_data_m) %in% "GSM907856"), ]

hip_pData   <- hip_pData[!(hip_pData$Source.Name %in% "GSM907837 1"),]
hip_pData   <- hip_pData[!(hip_pData$Source.Name %in% "GSM907856 1"),]


