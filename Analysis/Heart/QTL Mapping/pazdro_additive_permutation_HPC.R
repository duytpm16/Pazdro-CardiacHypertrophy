#####################################################################################################################################
#   This script is used to find significance threshold by permutation for each heart phenotype data sent by Robert Pazdro
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
col    <- args[1]
n_perm <- args[2]
cores  <- args[3]

load('pazdro_heart_phenotype_viewer_v1.Rdata')
dataset <- 'dataset.heart.phenotype'









### Extract log-transformed expression data
expr <- get(dataset)$data$norm


### Covar contains both sex and batch. Will use this for GDF11 and MSTN mapping.
#     Other phenotypes will be adjusted by sex (sex.covar) only.
covar <- get(dataset)$covar.matrix
sex.covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]













### Run permutation for phenotype
pheno_name <- colnames(expr)[i]


if(pheno_name %in% c('GDF11', 'MSTN')){
   perm <- scan1perm(genoprobs = genoprobs,
                     pheno     = expr[, pheno_name, drop = FALSE],
                     kinship   = K,
                     covar     = covar,
                     n_perm    = n_perm,
                     cores     = cores)
   
}else{
   perm <- scan1perm(genoprobs = genoprobs,
                     pheno     = expr[, pheno_name, drop = FALSE],
                     kinship   = K,
                     covar     = sex.covar,
                     n_perm    = n_perm,
                     cores     = cores)
}










### Save result
saveRDS(perm, file = paste0('pazdro_', pheno_name,'_permutation_', n_perm, '.rds'))

