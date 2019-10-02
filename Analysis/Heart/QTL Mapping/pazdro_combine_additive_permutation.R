### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2)




### Load data
load('~/Desktop/pazdro_heart_phenotype_viewer_v1.Rdata')




### Read in permuatations
body.weight    <- readRDS('pazdro_body.weight_additive_permutation_10000.rds')
heart.weight   <- readRDS('pazdro_heart.weight_additive_permutation_10000.rds')
tibia.length   <- readRDS('pazdro_tibia.length_additive_permutation_10000.rds')
wall.thickness <- readRDS('pazdro_wall.thickness_additive_permutation_10000.rds')
hw.bw.percent  <- readRDS('pazdro_hw.bw.percent_additive_permutation_10000.rds')
hw.bw.ratio    <- readRDS('pazdro_hw.bw.ratio_additive_permutation_10000.rds')
hw.tl.ratio    <- readRDS('pazdro_hw.tl.ratio_additive_permutation_10000.rds')
hw.reg.bw      <- readRDS('pazdro_hw.reg.bw_additive_permutation_10000.rds')
hw.reg.tl      <- readRDS('pazdro_hw.reg.tl_additive_permutation_10000.rds')
cc.area        <- readRDS('pazdro_cardiomyocyte.cross.sectional.area_additive_permutation_10000.rds')
perc.fibrosis  <- readRDS('pazdro_percent.fibrosis_additive_permutation_10000.rds')
gdf11          <- readRDS('pazdro_GDF11_additive_permutation_10000.rds')
mstn           <- readRDS('pazdro_MSTN_additive_permutation_10000.rds')
hw.adj.bw      <- readRDS('pazdro_hw.adj.bw_additive_permutation_10000.rds')
hw.adj.tl      <- readRDS('pazdro_hw.adj.tl_additive_permutation_10000.rds')








### Find significance at different alpha
body.weight.alpha    <- summary_scan1perm(object = body.weight,  alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
heart.weight.alpha   <- summary_scan1perm(object = heart.weight, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
tibia.length.alpha   <- summary_scan1perm(object = tibia.length, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
wall.thickness.alpha <- summary_scan1perm(object = wall.thickness, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.bw.percent.alpha  <- summary_scan1perm(object = hw.bw.percent, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.bw.ratio.alpha    <- summary_scan1perm(object = hw.bw.ratio, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.tl.ratio.alpha    <- summary_scan1perm(object = hw.tl.ratio, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.reg.bw.alpha      <- summary_scan1perm(object = hw.reg.bw, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.reg.tl.alpha      <- summary_scan1perm(object = hw.reg.tl, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
cc.area.alpha        <- summary_scan1perm(object = cc.area, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
perc.fibrosis.alpha  <- summary_scan1perm(object = perc.fibrosis, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
gdf11.alpha          <- summary_scan1perm(object = gdf11, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
mstn.alpha           <- summary_scan1perm(object = mstn, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.adj.bw.alpha      <- summary_scan1perm(object = hw.adj.bw, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))
hw.adj.tl.alpha      <- summary_scan1perm(object = hw.adj.tl, alpha = c(0.01, 0.05, 0.1, 0.2, 0.67))







### Combine
alpha <- cbind(body.weight.alpha, heart.weight.alpha, tibia.length.alpha,
               hw.bw.percent.alpha, hw.bw.ratio.alpha, hw.tl.ratio.alpha,
               gdf11.alpha, mstn.alpha, wall.thickness.alpha, cc.area.alpha, 
               perc.fibrosis.alpha, hw.adj.bw.alpha, hw.adj.tl.alpha,hw.reg.bw.alpha, 
               hw.reg.tl.alpha)
stopifnot(ncol(alpha) == 15)






dataset.heart.phenotype$perms$additive <- alpha




rm(list = ls()[-grep('dataset[.]|genoprobs|K|map|markers', ls())])
save.image('~/Desktop/pazdro_heart_phenotype_viewer_v1.Rdata')
