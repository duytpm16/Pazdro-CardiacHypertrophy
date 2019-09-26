#####################################################################################################################################
#   This script is used to run a sex-interaction scan on the heart phenotype data sent by Robert Pazdro
#
#
#
#   Input:
#     1.) QTL viewer formatted .Rdata file generated using 'pazdro_heart_phenotype_normalization.R'
#     2.) Additive scan LOD matrix
#
#
#   Output:
#     1.) .Rdata with tables in lod.peaks 'sex_int' index of dataset.heart.phenotype
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
dataset   <- 'dataset.heart.phenotype'
add_scan1 <- readRDS('pazdro_heart_phenotype_additive_scan1_matrix.rds')







### Extract log-transformed expression data
expr <- get(dataset)$data$norm
phys <- expr %>% as.data.frame() %>% select(-GDF11, -MSTN)
prot <- expr %>% as.data.frame() %>% select(GDF11, MSTN)


### Covar contains both sex and batch. Will use this for GDF11 and MSTN mapping.
#     Other phenotypes will be adjusted by sex (sex.covar) only.
covar <- get(dataset)$covar.matrix
sex.covar <- get(dataset)$covar.matrix[,'sex', drop = FALSE]










### Map GDF11 and MSTN proteins 
prot_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = prot,
                    kinship   = K,
                    addcovar  = covar,
                    intcovar  = sex.covar,
                    cores     = 5)



### Map heart physiological phenotypes
phys_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = phys,
                    kinship   = K,
                    addcovar  = sex.covar,
                    intcovar  = sex.covar,
                    cores     = 5)

### Combine scan1 outputs
scan1_out <- cbind(prot_scan1, phys_scan1)

stopifnot(colnames(add_scan1) == colnames(scan1_out))
diff <- scan1_out - add_scan1










### Find QTLs with LOD > 6
peaks <- find_peaks(scan1_output = diff, map = map, threshold = 6, drop = 1.5)
peaks <- peaks %>% 
          rename(data.name = lodcolumn, ci.lo = ci_lo, ci.hi = ci_hi) %>%
          left_join(x = ., y = markers[,c('marker.id', 'chr', 'pos')], by = c('chr', 'pos')) %>%
          select(data.name, marker.id, chr, pos, lod, ci.lo, ci.hi)










### Save to lod.peaks index as 'sex_int'
dataset.heart.phenotype$lod.peaks$sex_int <- peaks


### Save scan1 matrix
save(scan1_out, file = 'pazdro_heart_phenotype_sex_int_scan1_matrix.rds')
save(diff, file = 'pazdro_heart_phenotype_sex_effect_scan1_matrix.rds')





### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata') 
