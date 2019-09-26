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
expr <- get(dataset)$data$norm
phys <- expr %>% as.data.frame() %>% select(-GDF11, -MSTN)
prot <- expr %>% as.data.frame() %>% select(GDF11, MSTN)


### Covar contains both sex and batch. Will use this for GDF11 and MSTN mapping.
#     Other phenotypes will be adjusted by sex (sex.covar) only.
covar <- get(dataset)$covar.matrix
sex.covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]















### Estimate heritability per chromosome for each phenotype
chr_herit <- matrix(0, nrow = length(c(1:19,'X')), ncol = ncol(expr),
                    dimnames = list(paste0('chr', c(1:19,'X')), colnames(expr)))

for(i in colnames(chr_herit)){
    for(j in c(1:19,'X')){
        if(i %in% c('GDF11', 'MSTN')){
           chr_herit[paste0('chr', j), i] <- est_herit(pheno = prot[, i, drop = FALSE], kinship  = K[[j]], addcovar = covar)
           
        }else{
           chr_herit[paste0('chr', j), i] <- est_herit(pheno = phys[, i, drop = FALSE], kinship  = K[[j]], addcovar = sex.covar)           
        }
    }
}










### Overall heritability for each phenotype
K_overall <- calc_kinship(probs = genoprobs, type = 'overall', cores = 5)
overall_herit <- matrix(0, nrow = 1, ncol = ncol(expr),
                        dimnames = list('heritability', colnames(expr)))


for(i in colnames(overall_herit)){
    if(i %in% c('GDF11', 'MSTN')){
       overall_herit[, i] <- est_herit(pheno = prot[, i, drop = FALSE], kinship  = K_overall, addcovar = covar)
      
    }else{
       overall_herit[, i] <- est_herit(pheno = phys[, i, drop = FALSE], kinship  = K_overall, addcovar = sex.covar)           
    }
  }
}














### Store results to 'heritability index'
dataset.heart.phenotype$heritability$overall <- as.data.frame(overall_herit)
dataset.heart.phenotype$heritability$chr <- as.data.frame(chr_herit)









### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata') 
