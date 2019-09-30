#####################################################################################################################################
#   This script is used to find additive significance threshold by permutation for each heart phenotype data sent by Robert Pazdro
#       This runs in parallel on HPC
#
#
#   Input:
#     1.) QTL viewer formatted .Rdata file generated using 'pazdro_heart_phenotype_normalization.R'
#     2.) col    - which column of the expression matrix to grab
#     3.) n_perm - Number of permutations to run
#     4.) cores  - Number of cores to use
#
#   Output:
#     1.) .rds file containing the matrix (n_perm x 1) with maximum LOD score for each permutation run 
#
#
#
#   Author: Duy Pham
#   E-mail: duy.pham@jax.org
#####################################################################################################################################

### Options and libraries and seed (for reproducibility?)
options(stringsAsFactors = FALSE)
library(tidyverse)
library(qtl2)
set.seed(12345)







### Load data
args <- commandArgs(trailingOnly = TRUE)
col    <- as.numeric(args[1])
n_perm <- as.numeric(args[2])
cores  <- as.numeric(args[3])

load('pazdro_heart_phenotype_viewer_v1.Rdata')
dataset <- 'dataset.heart.phenotype'







### Run permutation for phenotype
pheno_name <- colnames(expr)[col]




### Extract log-transformed expression data
expr <- get(dataset)$data$norm


### Get covariates
if(pheno_name %in% c('GDF11', 'MSTN'){
   covar <- get(dataset)$covar.matrix[,c('sex','batch')]
}
if(pheno_name %in% c('hw.adj.bw'){
   covar <- get(dataset)$covar.matrix[,c('sex','body.weight')]
}
if(pheno_name %in% c('hw.adj.tl'){
   covar <- get(dataset)$covar.matrix[,c('sex','tibia.length')]
}
if(!pheno_name %in% c('GDF11', 'MSTN', 'hw.adj.bw', 'hw.adj.tl')){
   covar <-  get(dataset)$covar.matrix[,'sex', drop = FALSE]
}















perm <- scan1perm(genoprobs = genoprobs,
                  pheno     = expr[, pheno_name, drop = FALSE],
                  kinship   = K,
                  covar     = covar,
                  n_perm    = n_perm,
                  cores     = cores)











### Save result
saveRDS(perm, file = paste0('pazdro_', pheno_name,'_additive_permutation_', n_perm, '.rds'))

