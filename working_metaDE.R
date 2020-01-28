#Following Code implements Adaptively Weighted Fisher's method with One Sided Correct to join several datasets by using their p values. 

rm(list = ls())
head(pheno)
BiocManager::install("RSQLite")
library(annotate)lookUp
library(AnnotationDbi)
library(org.Hs.eg.db)
library(WGCNA)
library(AWFisher)
#load("male_meta_analysis.RData")

fisher_test_func <- function (p.values) 
{
  if (NCOL(p.values) == 1) 
    p.values = t(p.values)
  n = NCOL(p.values)
  if (n < 2 | n > 100) {
    stop("number of studies K has to be >= 2 and <= 100.")
  }
  out <- .C("AWpvalue", best_stat = rep(0, NROW(p.values)), 
            sum_weight = as.integer(rep(0, NROW(p.values))),
            weights = as.integer(rep(0, length(p.values))), 
            pval = t(p.values), nrow = as.integer(NROW(p.values)),ncol = as.integer(NCOL(p.values)))
  
  bestStat <- out$best_stat
  weights <- t(matrix(out$weights, ncol = NROW(p.values)))
  list(pvalues = aw.fisher.stat(bestStat, n), weights = weights)
}


# Lines 32-63 import individual datasets for our meta analysis
rootdir <- "C:/Users/axv851/Documents/Meta_Analysis_Pre_Processing/Analysis/male_pvalues"
setwd(rootdir)

load("male_meta_strat.RData")

eset_48350   <- read.table("E-GEOD-48350_eset.txt", sep = "\t")
#eset_5281    <- read.table("E-5281-metaDE.txt", sep = "\t")
eset_3157225 <- read.table("syn3157225_eset_final.txt")
#eset_29378   <- read.table("E-GEOD-29378_eset_final.txt", sep = "\t")
#eset_28146   <- read.table("E-GEOD-28146_eset_final.txt", sep = "\t")
#eset_36980   <- read.table("E_36980_eset.txt")


#colnames(eset_48350)   <- gsub("X","", colnames(eset_48350))
#colnames(eset_5281)    <- gsub("X","", colnames(eset_5281))
#colnames(eset_3157225) <- gsub("X","", colnames(eset_3157225))
#colnames(eset_28146)   <- gsub("X","", colnames(eset_28146))
#colnames(eset_29378)   <- gsub("X","", colnames(eset_29378))
#colnames(eset_36980)   <- gsub("X","", colnames(eset_36980))


#pheno_5281    <- eset_5281[,5316:5317]
#pheno_48350   <- read.csv("pheno_E48350.csv")
#pheno_3157225 <- eset_3157225[,2958:2959] 
#pheno_29378   <- read.csv("pdata_29378.csv")
#pheno_28146   <- read.csv("pheno_28146.csv")
#pheno_36980   <- eset_36980[ , 2144:2145 ]

pvals_3157225 <- read.csv("male_syn3157225.csv")
pvals_48350   <- read.csv("male_E-GEOD-48350.csv")
pvals_5281    <- read.csv("maleE-GEOD-5281.csv")
pvals_29378   <- read.csv("male_E-GEOD-29378.csv")
pvals_28146   <- read.csv("male_E-GEOD-28146.csv")
pvals_36980   <- read.csv("male_E-GEOD-36980.csv")


# We create lists of the individual results as well as create a matrix of weights, p values and standard errors to create
# a final table of logfc, logfc confidence intervals and q-values 
list_pvals   <- list(pvals_28146,pvals_29378,pvals_48350,pvals_5281,pvals_36980, pvals_3157225)
names(list_pvals) <- c("E28146", "E29378", "E48350", "E5281","E36980","syn3157225")



#pheno_48350    <- pheno_48350[,3:4]
#pheno_48350    <- pheno_48350[,c(2,1)]
#pheno_29378    <- pheno_29378[,3:4]
#pheno_28146    <- pheno_28146[,3:4]

#list_pheno_df  <- list(pheno_28146,pheno_29378,pheno_3157225,pheno_48350,pheno_5281,pheno_36980)
#names(list_pheno_df) <- c("pheno_28146","pheno_29378","pheno_3157225","pheno_48350","pheno_5281", "pheno_36980")
#names_c        <- c("case", "sex")



#list_pheno_df  <- lapply(list_pheno_df,setNames,names_c) 
#list_pheno_df$pheno_3157225$sex <- ifelse(list_pheno_df$pheno_3157225$sex == "F", "female", "male")
#list_pheno_m   <- lapply(list_pheno_df,function(x){subset(x,x$sex == "male")})

#eset_48350_f   <- eset_48350[rownames(eset_48350) %in% rownames(list_pheno_m$pheno_48350),]
#eset_3157725_f <- eset_3157225[rownames(eset_3157225) %in% rownames(list_pheno_m$pheno_syn3157225), ]
#eset_5281_f    <- eset_5281[rownames(eset_5281) %in% rownames(list_pheno_m$pheno_5281), ]
#eset_29378_f   <- subset(eset_29378, eset_29378$sex == "male")
#eset_36980_f   <- eset_36980[rownames(eset_36980) %in% rownames(list_pheno_m$pheno), ]
#eset_28146_f   <- eset_28146[rownames(eset_28146) %in% rownames(list_pheno_m$pheno_28146),]
#eset_28146_f   <- subset(eset_28146, eset_28146$sex == "male")
#n_col          <- ncol(eset_28146_f) - 2
#eset_28146_f   <- eset_28146_f[ , 1:n_col ]






list_genes        <- lapply(list_pvals, function(x){ as.character(lookUp(as.character(x$X),'org.Hs.eg','SYMBOL'))})
names(list_genes) <- names(list_pvals)

for (i in 1:6){
  list_pvals[[i]]$sym <- list_genes[[i]]
}

union_probes   <- multiUnion(list_genes)     #symbol_5281)

#common_probes  <- intersect(common_probes, symbol_5281)
common_probes  <- multiIntersect(list_genes)
#union_probes   <- union(union_probes, symbol_5281)

# pval_presence is a matrix representing the presence of a gene reading in the corresponding datasets where the rows correspond to the gene
# names and the columns correspond to the study, the weight assigned to studies with pval_presence value equal to false will be zero in Adaptively Weighted Fisher's 
pval_presence  <- data.frame(matrix(nrow = length(union_probes), ncol = 7))
colnames(pval_presence) <- c("gene", names(list_genes))
pval_presence$gene          <- union_probes
pval_presence$E48350        <- pval_presence$gene %in% list_genes$E48350
pval_presence$syn3157225    <- pval_presence$gene %in% list_genes$syn3157225
pval_presence$E28146        <- pval_presence$gene %in% list_genes$E28146
pval_presence$E5281         <- pval_presence$gene %in% list_genes$E5281
pval_presence$E29378        <- pval_presence$gene %in% list_genes$E29378
pval_presence$E36980        <- pval_presence$gene %in% list_genes$E36980

pval_presence_all   <- subset(pval_presence, pval_presence$all_studies == 6)


pval_mat           <- data.frame(matrix(nrow = length(union_probes), ncol = 7))
colnames(pval_mat) <- colnames(pval_presence)

# Pval_mat is a dataframe with the rows representing genes and the columns representing studies 
# If a gene "X" is missing from study # "N", we insert the P value 1 as the weight for studies with pvalue = 1 is zero  
# in the calculation of the combined P value
pval_mat$gene      <- pval_presence$gene

pval_mat$syn3157225   <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$syn3157225 == TRUE, list_pvals$syn3157225[list_pvals$syn3157225$sym == x,]$adj.P.Val,1)})
pval_mat$E48350       <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$E48350     == TRUE, list_pvals$E48350[list_pvals$E48350$sym == x,]$adj.P.Val,1)})
pval_mat$E5281        <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$E5281      == TRUE, list_pvals$E5281[list_pvals$E5281$sym == x,]$adj.P.Val,1)})
pval_mat$E28146       <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$E28146     == TRUE, list_pvals$E28146[list_pvals$E28146$sym == x,]$adj.P.Val,1)})
pval_mat$E29378       <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$E29378     == TRUE, list_pvals$E29378[list_pvals$E29378$sym == x,]$adj.P.Val,1)})
pval_mat$E36980       <- lapply(pval_presence$gene, function(x){ ifelse(pval_presence[pval_presence$gene == x,]$E36980     == TRUE, list_pvals$E36980[list_pvals$E36980$sym == x,]$adj.P.Val,1)})


pval_mat <- pval_mat[,8:ncol(pval_mat)]

pval_mat_all <- subset(pval_mat, pval_mat$gene %in% pval_presence_all$gene)





pval_mat_int <- as.data.frame(pval_mat[,2:7], ncol = 6)
pval_mat_int_x <- as.matrix(pval_mat_int)
pval_mat_int_x <- apply(pval_mat_int_x,2,as.numeric)
rownames(pval_mat_int) <- pval_mat[,1]
rownames(pval_mat_int_x) <- pval_mat$gene
#save.image("male_meta_strat.RData")


meta_test <- AWFisher_pvalue(pval_mat_int_x)

wt_table           <- meta_test$weights 
colnames(wt_table) <- colnames(pval_mat_int)
rownames(wt_table) <- rownames(pval_mat_int)
wt_table           <- data.frame(wt_table)
qvalue             <- p.adjust(meta_test$pvalues, method = "BH")

qvalue_df <- as.data.frame(cbind(pval_mat$gene, as.character(qvalue)))
colnames(qvalue_df) <- c("gene", "qvalue")
#qvalue_df$qvalue <- as.double(qvalue_df$qvalue)



wt_table$gene         <- pval_mat$gene

# After extracting the weights of the studies for individual probes and the p values, we need to get standard errors and the the logfc values to 
# to get results from final meta.summary which combines results using weights

# The following lines of code populate the Log Fold Change matrix based on the weights table("wt_table")
logfc_mat              <- data.frame(matrix(nrow = length(union_probes), ncol = 7 ))
colnames(logfc_mat)    <- colnames(pval_mat)
rownames(logfc_mat)    <- rownames(pval_mat)

logfc_mat$gene            <- wt_table$gene
logfc_mat$syn3157225      <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$syn3157225 == 1, list_pvals$syn3157225[list_pvals$syn3157225$sym == x,]$logFC ,1)})
logfc_mat$E48350          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E48350     == 1, list_pvals$E48350[list_pvals$E48350$sym == x,]$logFC,1)})
logfc_mat$E5281           <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E5281      == 1, list_pvals$E5281[list_pvals$E5281$sym == x,]$logFC,1)})
logfc_mat$E28146          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E28146     == 1, list_pvals$E28146[list_pvals$E28146$sym == x,]$logFC,1)})
logfc_mat$E29378          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E29378     == 1, list_pvals$E29378[list_pvals$E29378$sym == x,]$logFC,1)})
logfc_mat$E36980          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene ==x,]$E36980 == 1,      list_pvals$E36980[list_pvals$E36980$sym == x, ]$logFC,1)})


# The following lines of code populate the Standard Error  matrix based on the weights table("wt_table")
se_mat     <- data.frame(matrix(nrow = length(union_probes), ncol = 7))
colnames(se_mat)  <- colnames(pval_mat)
rownames(se_mat)  <- rownames(pval_mat)

se_mat$gene    <- wt_table$gene

se_mat$syn3157225      <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$syn3157225 == 1,    (list_pvals$syn3157225[list_pvals$syn3157225$sym == x,]$CI.R - list_pvals$syn3157225[list_pvals$syn3157225$sym == x,]$CI.L)*0.5 ,1)})
se_mat$E48350          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E48350 == 1,        (list_pvals$E48350[list_pvals$E48350$sym == x,]$CI.R - list_pvals$E48350[list_pvals$E48350$sym == x,]$CI.L)*0.5 ,1)})
se_mat$E5281           <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E5281 == 1,         (list_pvals$E5281[list_pvals$E5281$sym == x,]$CI.R - list_pvals$E5281[list_pvals$E5281$sym == x,]$CI.L)*0.5 ,1)})
se_mat$E28146          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E28146 == 1,        (list_pvals$E28146[list_pvals$E28146$sym == x,]$CI.R - list_pvals$E28146[list_pvals$E28146$sym == x,]$CI.L)*0.5 ,1)})
se_mat$E29378          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E29378 == 1,        (list_pvals$E29378[list_pvals$E29378$sym == x,]$CI.R - list_pvals$E29378[list_pvals$E29378$sym == x,]$CI.L)*0.5 ,1)})
se_mat$E36980          <- lapply(wt_table$gene, function(x){ ifelse(wt_table[wt_table$gene == x,]$E36980 == 1,        (list_pvals$E36980[list_pvals$E36980$sym == x,]$CI.R - list_pvals$E36980[list_pvals$E36980$sym == x,]$CI.L)*0.5 ,1)})

se_mat_all    <-  subset(se_mat,se_mat$gene %in% pval_presence_all$gene)
logfc_mat_all <-  subset(logfc_mat, logfc_mat$gene %in% pval_presence_all$gene)

meta_gen_method  <- metagen(TE = logfc_mat_all[1,2:ncol(logfc_mat_all)], seTE = se_mat_all[1,2:ncol(se_mat_all)])

meta_sum <- lapply(wt_table$gene, function(x){ meta.summaries(d = as.vector(as.numeric(unlist(logfc_mat[logfc_mat$gene == x,2:7]))),
                                                              se = as.vector(as.numeric(unlist(se_mat[se_mat$gene ==x,2:7]))), method = "fixed", 
                                                              weights = as.vector(as.numeric(unlist(wt_table[wt_table$gene == x,1:6]))))})   # Compute meta_sum for all the genes in our gene set

names(meta_sum)  <- wt_table$gene
meta_summarry    <- lapply(meta_sum, function(x){summary(x)})
meta_rest        <- lapply(meta_sum, function(x){summary(x)$summci})
 
meta_rest_df           <- do.call("rbind", meta_rest)
colnames(meta_rest_df) <- c("Lower_CI","log_fc","Upper_CI")
meta_rest_df           <- as.data.frame(cbind(meta_rest_df,as.numeric(meta_test$pvalues)))
meta_rest_df           <- meta_rest_df[,c(2,1,3,4)]
colnames(meta_rest_df) <- c("log_fc","lower_ci","upper_ci", "adjusted_pval_with_AW")
meta_rest_df    <- as.data.frame(meta_rest_df)

write.csv(meta_rest_df, "female_final_meta_analysis.csv")



