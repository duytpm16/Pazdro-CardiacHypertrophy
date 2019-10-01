#####################################################################################################################################
#   This script is used run additive scan on the heart phenotype data sent by Robert Pazdro
#
#
#
#   Input:
#     1.) QTL viewer formatted .Rdata file generated using 'pazdro_heart_phenotype_normalization.R'
#
#
#   Output:
#     1.) .Rdata with tables in lod.peaks 'additive' index of dataset.heart.phenotype
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
phys <- expr %>% as.data.frame() %>% select(-GDF11, -MSTN, -hw.adj.bw, -hw.adj.tl)
prot <- expr %>% as.data.frame() %>% select(GDF11, MSTN)
hw.adj.bw <- expr %>% as.data.frame() %>% select(hw.adj.bw)
hw.adj.tl <- expr %>% as.data.frame() %>% select(hw.adj.tl)

### Get covariates
prot.covar <- get(dataset)$covar.matrix[,c('sex','batch'), drop = FALSE]
sex.covar  <- get(dataset)$covar.matrix[,'sex', drop = FALSE]
bw.covar   <- get(dataset)$covar.matrix[,c('sex','body.weight'), drop = FALSE]
tl.covar   <- get(dataset)$covar.matrix[,c('sex','tibia.length'), drop = FALSE]












### Map GDF11 and MSTN proteins 
prot_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = prot,
                    kinship   = K,
                    addcovar  = prot.covar,
                    cores     = 5)



### Map heart physiological phenotypes
phys_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = phys,
                    kinship   = K,
                    addcovar  = sex.covar,
                    cores     = 5)

### Map heart weight adjusting for bodyweight
hw.bw_scan1 <- scan1(genoprobs = genoprobs,
                     pheno     = hw.adj.bw,
                     kinship   = K,
                     addcovar  = bw.covar,
                     cores     = 5)


### Map heart weight adjusting for bodyweight
hw.tl_scan1 <- scan1(genoprobs = genoprobs,
                     pheno     = hw.adj.tl,
                     kinship   = K,
                     addcovar  = tl.covar,
                     cores     = 5)

### Combine scan1 outputs
scan1_out <- cbind(prot_scan1, phys_scan1, hw.bw_scan1, hw.tl_scan1)












### Find QTLs with LOD > 5.5
peaks <- find_peaks(scan1_output = scan1_out, map = map, threshold = 5.5, drop = 1.5)
peaks <- peaks %>% 
           rename(data.name = lodcolumn, ci.lo = ci_lo, ci.hi = ci_hi) %>%
           left_join(x = ., y = markers[,c('marker.id', 'chr', 'pos')], by = c('chr', 'pos')) %>%
           select(data.name, marker.id, chr, pos, lod, ci.lo, ci.hi)










### Add allele effects
peaks[,LETTERS[1:8]] <- 0
for(i in 1:nrow(peaks)){
  
    gp <- genoprobs[,peaks$chr[i]]
    gp[[1]] <- gp[[1]][,,peaks$marker.id[i], drop = FALSE]
  
    covar <- switch (peaks$data.name[i],
                      'GDF11'= prot.covar,
                      'MSTN' = prot.covar,
                      'hw.adj.bw' = bw.covar,
                      'hw.adj.tl' = tl.covar,
                      sex.covar)
 
    peaks[i,LETTERS[1:8]] <- scan1blup(genoprobs = gp,
                                       pheno     = expr[,peaks$data.name[i], drop = FALSE],
                                       kinship   = K[[peaks$chr[i]]],
                                       addcovar  = covar)[,LETTERS[1:8]]      
    
}












### Save to lod.peaks index as 'additive'
dataset.heart.phenotype$lod.peaks$additive <- peaks


### Save scan1 matrix
saveRDS(scan1_out, file = 'pazdro_heart_phenotype_additive_scan1_matrix.rds')






### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata') 
