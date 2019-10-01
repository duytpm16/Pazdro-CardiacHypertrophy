#####################################################################################################################################
#   This script is used to estimate heritability on the heart phenotype data sent by Robert Pazdro
#
#
#
#   Input:
#     1.) QTL viewer formatted .Rdata file generated using 'pazdro_heart_phenotype_normalization.R'
#
#
#   Output:
#     1.) .Rdata with 'herit' list in dataset.heart.phenotype list
#           Contains: 'overall' which is overall heritability estimates
#                     'chr' which is heritability estimates for each chromosome
#
#
#
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#####################################################################################################################################

### Options and libraries
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)








### Load data
load('pazdro_heart_phenotype_viewer_v1.Rdata')
dataset <- 'dataset.heart.phenotype'








### Extract log-transformed expression data
expr <- get(dataset)$data$log

### Get covariates
prot.covar <- get(dataset)$covar.matrix[,c('sex','batch'), drop = FALSE]
sex.covar  <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
bw.covar   <- get(dataset)$covar.matrix[,c('sex','body.weight'), drop = FALSE]
tl.covar   <- get(dataset)$covar.matrix[,c('sex','tibia.length'), drop = FALSE]















### Estimate heritability per chromosome for each phenotype
chr_herit <- matrix(0, nrow = length(c(1:19,'X')), ncol = ncol(expr),
                    dimnames = list(paste0('chr', c(1:19,'X')), colnames(expr)))

for(i in colnames(chr_herit)){
    for(j in c(1:19,'X')){
    
        covar <- switch(i,
                        'GDF11'= prot.covar,
                        'MSTN' = prot.covar,
                        'hw.adj.bw' = bw.covar,
                        'hw.adj.tl' = tl.covar,
                        sex.covar)
    
        chr_herit[paste0('chr', j), i] <- est_herit(pheno = expr[, i, drop = FALSE], kinship  = K[[j]], addcovar = covar)           
    }
}










### Overall heritability for each phenotype
K_overall <- calc_kinship(probs = genoprobs, type = 'overall', cores = 5)
overall_herit <- matrix(0, nrow = 1, ncol = ncol(expr),
                        dimnames = list('heritability', colnames(expr)))


for(i in colnames(overall_herit)){

    covar <- switch (i,
                     'GDF11'= prot.covar,
                     'MSTN' = prot.covar,
                     'hw.adj.bw' = bw.covar,
                     'hw.adj.tl' = tl.covar,
                     sex.covar)
  
    overall_herit[, i] <- est_herit(pheno = expr[, i, drop = FALSE], kinship  = K_overall, addcovar = covar)           
}














### Store results to 'heritability index'
dataset.heart.phenotype$heritability$overall <- as.data.frame(overall_herit)
dataset.heart.phenotype$heritability$chr <- as.data.frame(chr_herit)









### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.RData') 
