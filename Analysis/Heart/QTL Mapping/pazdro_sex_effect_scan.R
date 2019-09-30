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
expr <- get(dataset)$data$log
phys <- expr %>% as.data.frame() %>% select(-GDF11, -MSTN, -hw.adj.bw, -hw.adj.tl)
prot <- expr %>% as.data.frame() %>% select(GDF11, MSTN)
hw.adj.bw <- expr %>% as.data.frame() %>% select(hw.adj.bw)
hw.adj.tl <- expr %>% as.data.frame() %>% select(hw.adj.tl)

### Covar contains both sex and batch. Will use this for GDF11 and MSTN mapping.
#     Other phenotypes will be adjusted by sex (sex.covar) only.
prot.covar <- get(dataset)$covar.matrix[,c('sex','batch'), drop = FALSE]
sex.covar  <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
bw.covar   <- get(dataset)$covar.matrix[,c('sex','body.weight'), drop = FALSE]
tl.covar   <- get(dataset)$covar.matrix[,c('sex','tibia.length'), drop = FALSE]










### Map GDF11 and MSTN proteins 
prot_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = prot,
                    kinship   = K,
                    addcovar  = prot.covar,
                    intcovar  = sex.covar,
                    cores     = 5)

### Map heart physiological phenotypes
phys_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = phys,
                    kinship   = K,
                    addcovar  = sex.covar,
                    intcovar  = sex.covar,
                    cores     = 5)

### Map heart weight adjusting for bodyweight
hw.bw_scan1 <- scan1(genoprobs = genoprobs,
                     pheno     = hw.adj.bw,
                     kinship   = K,
                     addcovar  = bw.covar,
                     intcovar  = sex.covar,
                     cores     = 5)


### Map heart weight adjusting for bodyweight
hw.tl_scan1 <- scan1(genoprobs = genoprobs,
                     pheno     = hw.adj.tl,
                     kinship   = K,
                     addcovar  = tl.covar,
                     intcovar  = sex.covar,
                     cores     = 5)

### Combine scan1 outputs
scan1_out <- cbind(prot_scan1, phys_scan1, hw.bw_scan1, hw.tl_scan1)







### Get sex-effects QTLs
stopifnot(colnames(add_scan1) == colnames(scan1_out))
stopifnot(rownames(add_scan1) == rownames(scan1_out))
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
saveRDS(scan1_out, file = 'pazdro_heart_phenotype_sex_int_scan1_matrix.rds')
saveRDS(diff, file = 'pazdro_heart_phenotype_sex_effect_scan1_matrix.rds')





### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata') 
