#####################################################################################################################################
#   This script is used to normalize the heart phenotype data sent by Robert Pazdro
#     Data are then formatted for QTL viewer
#
#
#   Genoprobs, map, markers, K were obtained on cadillac HPC at /projects/churchill-lab/data/Pazdro/genotypes/genoprobs/gigaMUGA/qtl2/
#
#
#
#   Input:
#     1.) .Rdata file on line 6 above
#     2.) phenotype data located on cadillac - /projects/churchill-lab/data/Pazdro/phenotypes/heart/
#
#
#   Output:
#     1.) .Rdata in QTL viewer format
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
#    1.) GigaMUGA genoprobs with 137k markers: /projects/churchill-lab/data/Pazdro/genotypes/genoprobs/gigaMUGA/qtl2/
#    2.) Pazdro's data: /projects/churchill-lab/data/Pazdro/phenotypes/heart/
load('~/Desktop/Pazdro Cardiac Hypertrophy/Genotypes/UNC/pazdro_gigaMUGA_genoprobs_qtl2.Rdata')
data <- readRDS('~/Desktop/Pazdro Cardiac Hypertrophy/Phenotype/Modified Phenotype/pazdro_heart_phenotype_cleaned_spreadsheet_plus_echo_2_14_2019.rds')









### Create samples dataframe
annot.sample <- data %>% 
                  select(Mouse.Number, Sex, DOB, Tissue.Collection.Date, `Batch.GDF11/MSTN`, Body.Weight.g, Tibia.Length.mm) %>%
                  dplyr::rename(mouse.id = Mouse.Number,
                                batch = `Batch.GDF11/MSTN`,
                                body.weight = Body.Weight.g,
                                tibia.length = Tibia.Length.mm) %>%
                  mutate(Sex = factor(Sex), batch = factor(batch),
                         body.weight = log(body.weight), tibia.length = log(tibia.length)) %>%
                  `colnames<-`(tolower(colnames(.))) %>%
                  select(mouse.id, sex, batch, dob, tissue.collection.date)










### Raw data formatting
data <- data %>%
          select(-Sex, -DOB, -Tissue.Collection.Date, -`Batch.GDF11/MSTN`) %>%
          select(-grep('LV|Fraction|Stroke|Cardiac.Output', colnames(.))) %>%
          dplyr::rename(body.weight    = Body.Weight.g,
                        heart.weight   = Heart.Weight.mg,
                        tibia.length   = Tibia.Length.mm,
                        hw.bw.percent  = `HW/BW.%`,
                        hw.bw.ratio    = `HW/BW.mg/g`,
                        hw.tl.ratio    = `HW/TL.mg/mm`,
                        GDF11          = `GDF11.ng/mL`,
                        MSTN           = `MSTN.ng/mL`,
                        wall.thickness = `Wall.Thickness.µm`,
                        cardiomyocyte.cross.sectional.area = `Cardiomyocyte.Cross.Sectional.Area.µm2`,
                        percent.fibrosis = `Percent.Fibrosis.% area`) %>%
          mutate(hw.adj.bw = heart.weight,
                 hw.adj.tl = heart.weight) %>%
          remove_rownames() %>%
          column_to_rownames('Mouse.Number')








### Normalize data
norm <- log(data)


### Temporary data frame for regression
temp <- cbind(model.matrix(~sex, data = annot.sample)[,-1,drop = FALSE], norm)

### Regress out body weight in heart weight data
fit  <- lm(heart.weight ~ sexM + body.weight, data = temp, na.action = na.exclude)
coef <- coefficients(fit)
res  <- residuals(fit)
hw.reg.bw <- t(coef[1] + (coef[2] %*% temp[,names(coef)[2]]) + res)

### Regress out tibia length in heart weight data
fit  <- lm(heart.weight ~ sexM + tibia.length, data = temp, na.action = na.exclude)
coef <- coefficients(fit)
res  <- residuals(fit)
hw.reg.tl <- t(coef[1] + (coef[2] %*% temp[,names(coef)[2]]) + res)



### Combined adjusted values to normalized data
norm <- cbind(norm, hw.reg.bw, hw.reg.tl)






### RankZ
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
} # rankZ()


rz = apply(norm, 2, rankZ)














### Covariate matrix
covar.matrix <- model.matrix(~ sex + batch, data = annot.sample)[, -1, drop = FALSE]
colnames(covar.matrix) <- c('sex', 'batch')
rownames(covar.matrix) <- annot.sample$mouse.id

stopifnot(rownames(norm) == rownames(covar.matrix))
covar.matrix <- cbind(covar.matrix, norm[,c('body.weight','tibia.length')])







### Covar info
covar.info <- data.frame(sample.column = c('sex', 'batch', 'body.weight', 'tibia.length'),
                         covar.column  = c('sex', 'batch', 'body.weight', 'tibia.length'),
                         display.name  = c('Sex', 'Batch', 'Body Weight', 'Tibia Length'),
                         interactive   = c(TRUE, FALSE, FALSE, FALSE),
                         primary       = c(TRUE, FALSE, FALSE, FALSE),
                         lod.peaks     = c('sex_int', NA, NA, NA))
















### Creating annot.phenotype dataframe
sample.col <- colnames(annot.sample)[-c(1:2)]
annot.phenotype <- data.frame(data.name      = c(sample.col, colnames(norm)),
                              short.name     = c(sample.col, colnames(norm)),
                              R.name         = c(sample.col, colnames(norm)),
                              description    = c('Mouse identifier', 'Sex of mouse: Female (F) or Male (M)', 'Mass spectrometry batch group for GDF11 and MSTN', 'Mouse data-of-birth', 'Tissue collection date',
                                                 'Body weight measured in grams', 'Heart weight measured in milligrams', ' Tibia length measured in millimeters', 
                                                 '(Heart weight / 100) / (Body weight / 100)', 'Heart weight to body weight ratio',
                                                 'Heart weight to tibia length ratio', 'GDF11 levels', 'MSTN levels', 'Thickness of heart wall measured in micrometer',
                                                 'Cardiomyocyte cross sectional area', 'Fibrosis percentage', 'Heart weight conditioned on body weight',
                                                 'Heart weight conditioned on tibia length', 'Heart weight with body weight regressed out', 'Heart weight with tibia length regressed out'),
                              units          = c(rep(NA, length(sample.col)), 'g', 'mg', 'mm', '%', 'mg/g', 'mg/mm', 'ng/mL', 'ng/mL', 'µm', 'µm2', '% area', NA, NA, NA, NA),
                              category       = c(rep('Demographic', length(sample.col)), rep('Phenotype', ncol(norm))),
                              R.category     = c(rep('Demographic', length(sample.col)), rep('Phenotype', ncol(norm))),
                              is.id          = c(TRUE, rep(FALSE, length(sample.col)) + ncol(norm) - 1)),
                              is.numeric     = c(rep(FALSE, length(sample.col)), rep(TRUE, ncol(norm))),
                              is.date        = c(FALSE, FALSE , FALSE, TRUE, TRUE, rep(FALSE, ncol(norm))),
                              is.factor      = c(FALSE, TRUE, TRUE, FALSE , FALSE, rep(FALSE, ncol(norm))),
                              factor.levels  = c(NA, 'F:M', '1:2', NA, NA, rep(NA, ncol(norm))),
                              is.covar       = c(FALSE, TRUE, TRUE, FALSE , FALSE, rep(FALSE, ncol(norm))),
                              is.pheno       = c(rep(FALSE, length(sample.col)), rep(TRUE, ncol(norm))),
                              is.derived     = c(rep(FALSE, length(sample.col)), FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, rep(FALSE, ncol(norm) - 8), TRUE, TRUE),
                              omit           = FALSE,
                              use.covar      = c(rep(NA, length(sample.col)), rep('sex', 6), 'sex:batch', 'sex:batch', rep('sex', 3), 'sex:body.weight', 'sex:tibia.length', 'sex', 'sex'),
                              transformation = c(rep(NA, length(sample.col)), rep('log', ncol(norm))),
                              adj.pheno      = c(rep(NA, length(sample.col)), rep(NA, ncol(norm) - 4), 'body.weight', 'tibia.length', NA, NA))

















### Format for QTL viewer
dataset.heart.phenotype <- list(annot.phenotype = as_tibble(annot.phenotype),
                                annot.samples   = as_tibble(annot.sample),
                                covar.matrix    = as.matrix(covar.matrix),
                                covar.info      = as_tibble(covar.info),
                                data            = list(raw  = as.matrix(data),
                                                       log  = as.matrix(norm),
                                                       rz   = as.matrix(rz)),
                                datatype        = 'phenotype',
                                display.name    = 'Heart physiological phenotypes',
                                lod.peaks       = list())







### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.RData')          
