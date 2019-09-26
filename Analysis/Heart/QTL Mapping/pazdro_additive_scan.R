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
                    cores     = 5)



### Map heart physiological phenotypes
phys_scan1 <- scan1(genoprobs = genoprobs,
                    pheno     = phys,
                    kinship   = K,
                    addcovar  = sex.covar,
                    cores     = 5)

### Combine scan1 outputs
scan1_out <- cbind(prot_scan1, phys_scan1)












### Find QTLs with LOD > 5.5
peaks <- find_peaks(scan1_output = scan1_out, map = map, threshold = 5.5, drop = 1.5)
peaks <- peaks %>% 
           rename(data.name = lodcolumn, ci.lo = ci_lo, ci.hi = ci_hi) %>%
           left_join(x = ., y = markers[,c('marker.id', 'chr', 'pos')], by = c('chr', 'pos')) %>%
           select(data.name, marker.id, chr, pos, lod, ci.lo, ci.hi)
peaks[,LETTERS[1:8]] <- 0










### Add allele effects
for(i in 1:nrow(peaks)){
  
    gp <- genoprobs[,peaks$chr[i]]
    gp[[1]] <- gp[[1]][,,peaks$marker.id[i], drop = FALSE]
    
    if(peaks$data.name[i] %in% c('GDF11','MSTN')){
       peaks[i,LETTERS[1:8]] <- scan1blup(genoprobs = gp,
                                          pheno     = prot[,peaks$data.name[i], drop = FALSE],
                                          kinship   = K[[peaks$chr[i]]],
                                          addcovar  = covar)[,LETTERS[1:8]]
    }else{
       peaks[i,LETTERS[1:8]] <- scan1blup(genoprobs = gp,
                                          pheno     = phys[,peaks$data.name[i], drop = FALSE],
                                          kinship   = K[[peaks$chr[i]]],
                                          addcovar  = sex.covar)[,LETTERS[1:8]]      
      
    }
}






### Save to lod.peaks index as 'additive'
dataset.heart.phenotype$lod.peaks$additive <- peaks


### Save scan1 matrix
saveRDS(scan1_out, file = 'pazdro_heart_phenotype_additive_scan1_matrix.rds')






### Save .Rdata
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata') 
