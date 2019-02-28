### Options and libraries
options(stringsAsFactors = FALSE)
library(qtl2convert)
library(qtl2)
library(tidyverse)
setwd('~/Desktop')





### Load UNC markers annotation
### Load genoprobs for both wave
load("~/Desktop/Pazdro Cardiac Hypertrophy/Genotypes/UNC/snps.gigamuga.Rdata")
wave1 <- readRDS("~/Desktop/Pazdro Cardiac Hypertrophy/Genotypes/UNC/Robert_Pazdro_Wave1__genoprobs_8state.Rdata")
wave2 <- readRDS("~/Desktop/Pazdro Cardiac Hypertrophy/Genotypes/UNC/Robert_Pazdro_Wave2__genoprobs_8state.Rdata")



wave1 <- probs_qtl2_to_doqtl(wave1)
wave2 <- probs_qtl2_to_doqtl(wave2)













### Making a dataframe of 'DO-###' and original sample names for wave1 and wave2
wave1_name <- dimnames(wave1)[[1]] %>%
                                   data.frame(orig_genoprobs_name = .) %>%
                                   separate(orig_genoprobs_name, into = c('Size', 'File', 'DO', 'End'), sep = '_')  %>%
                                   separate(DO, into = c('DO', 'Number'), sep = '-') %>%
                                   dplyr::mutate(Number = as.numeric(Number)) %>%
                                   arrange(Number) %>%
                                   unite(col = 'DO', DO, Number, sep = '-') %>%
                                   group_by(DO) %>%
                                   dplyr::mutate(orig_genoprobs_name = paste0(c(Size, File, DO, End), collapse = '_')) %>%
                                   dplyr::select(-Size, -File, -End)


wave2_name <- dimnames(wave2)[[1]] %>%
                                   data.frame(orig_genoprobs_name = .) %>%
                                   separate(orig_genoprobs_name, into = c('Wave', 'DO', 'End'), sep = '_')  %>%
                                   separate(DO, into = c('DO', 'Number'), sep = '-') %>%
                                   dplyr::mutate(Number = as.numeric(Number)) %>%
                                   arrange(Number) %>%
                                   unite(col = 'DO', DO, Number, sep = '-') %>%
                                   group_by(DO) %>%
                                   dplyr::mutate(orig_genoprobs_name = paste0(c(Wave, DO, End), collapse = '_')) %>%
                                   dplyr::select(-Wave, -End)
wave_name  <- rbind(wave1_name, wave2_name)











### Create new array by getting genoprobs from wave1 and wave2 and arranging them by DO number
stopifnot(dim(wave1)[3] == dim(wave2)[3])
stopifnot(dimnames(wave1)[[3]] == dimnames(wave1)[[3]])
new_array <- array(0, dim = c(219, 8 , dim(wave1)[3]), dimnames = list(wave_name$DO, LETTERS[1:8], dimnames(wave1)[[3]]))

for(i in 1:dim(wave1)[1]){
    new_array[i,,] <- wave1[wave_name$orig_genoprobs_name[i],,]
}
for(i in 124:219){
    new_array[i,,] <- wave2[wave_name$orig_genoprobs_name[i],,]
}








### Create new marker map
markers <- snps %>%
                dplyr::mutate(chr = gsub('chr', '', chr, fixed = FALSE)) %>%
                filter(chr %in% c(1:19,'X')) %>%
                dplyr::mutate(bp  = as.numeric(pos),
                              chr = factor(chr, levels = c(1:19,'X')),
                              pos = as.numeric(pos) / 1e6) %>%
                select(chr, marker, cM, pos, rsID) %>%
                filter(marker %in% dimnames(new_array)[[3]])



### Create new marker map list
map <- map_df_to_list(map = markers,
                      chr_column = 'chr',
                      pos_column = 'pos')



### Convert new array to qtl2 genoprobs format
genoprobs <- probs_doqtl_to_qtl2(new_array, map = markers, pos_column = 'pos')




### Compute Kinship
K <- calc_kinship(probs = genoprobs,
                  type = 'loco',
                  cores = 10)


save(genoprobs, K, map , markers, file='pazdro_combined_UNC_genoprobs_with_map_K.RData')


