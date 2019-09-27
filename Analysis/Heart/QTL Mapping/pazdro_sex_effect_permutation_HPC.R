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
col    <- as.numeric(args[1])
n_perm <- as.numeric(args[2])
cores  <- as.numeric(args[3])

load('pazdro_heart_phenotype_viewer_v1.Rdata')
dataset <- 'dataset.heart.phenotype'







### Extract log-transformed expression data
expr <- get(dataset)$data$norm






### Get phenotype
pheno_name <- colnames(expr)[col]





### Get covariates
if(pheno_name %in% c('GDF11', 'MSTN')){
   covar <- get(dataset)$covar.matrix
   sex.covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
}else{
   covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
   sex.covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
}












### Generate permutations for genoprobs
orig_id <- dimnames(genoprobs[[1]])[[1]]

n_samp <- length(orig_id)
perms  <- lapply(1:n_perm, function(x) orig_id[sample(x = 1:n_samp, size = n_samp, replace = FALSE)])













### Additive permutation scans begin
add_scan1_output <- list()
for(i in 1:length(perms)){
  
    gp <- lapply(genoprobs, function(x){dimnames(x)[[1]] <- perms[[i]]; x})
    attributes(gp) <- attributes(genoprobs)
   
    add_scan1_output[[i]] <- scan1(genoprobs = gp, 
                                   pheno     = expr[, pheno_name,drop = FALSE],
                                   kinship   = K,
                                   addcovar  = covar,
                                   intcovar  = NULL,
                                   cores     = cores)
  
    print(paste0('additive ', i, ' of ', length(perms)))
}


### Transform to dataframe
add_scan1_output <- do.call(cbind, add_scan1_output)
colnames(add_scan1_output) <- paste0(colnames(add_scan1_output), '.perm.', 1:ncol(add_scan1_output))
saveRDS(object = add_scan1_output, file = paste0('pazdro_additive_manualPerm_LOD_matrix_', n_perm, '.rds'))





### Sex interaction permutation scans begin
int_scan1_output <- list()
for(i in 1:length(perms)){
  
    gp <- lapply(genoprobs, function(x){dimnames(x)[[1]] <- perms[[i]]; x})
    attributes(gp) <- attributes(genoprobs)
    
    int_scan1_output[[i]] <- scan1(genoprobs = gp, 
                                   pheno     = expr[, pheno_name,drop = FALSE],
                                   kinship   = K,
                                   addcovar  = covar,
                                   intcovar  = sex.covar,
                                   cores     = cores)
    
    print(paste0('interaction ', i, ' of ', length(perms)))    
}


### Transform to dataframe
int_scan1_output <- do.call(cbind, int_scan1_output)
colnames(int_scan1_output) <- paste0(colnames(int_scan1_output), '.perm.', 1:ncol(int_scan1_output))
saveRDS(object = int_scan1_output, file = paste0('pazdro_sex_int_manualPerm_LOD_matrix_', n_perm, '.rds'))














### Get sex-effect matrix
effect_output <- int_scan1_output - add_scan1_output
stopifnot(all(effect_output > 0))



### Find max for each permutation run
effect_max <- data.frame(apply(effect_output, 2, max))
colnames(effect_max) <- pheno_name










### Save result
saveRDS(effect_max, file = paste0('pazdro_', pheno_name,'_sex_effect_permutation_', n_perm, '.rds'))
