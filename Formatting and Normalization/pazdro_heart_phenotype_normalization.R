#####################################################################################################################################
#   This script is used normalize the heart phenotype data sent by Robert Pazdro
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
                  select(Mouse.Number, Sex, DOB, Tissue.Collection.Date, `Batch.GDF11/MSTN`) %>%
                  rename(mouse.id = Mouse.Number,
                         batch = `Batch.GDF11/MSTN`) %>%
                  mutate(Sex = factor(Sex), batch = factor(batch)) %>%
                  `colnames<-`(tolower(colnames(.))) %>%
                  select(mouse.id, sex, batch, dob, tissue.collection.date)










### Raw data formatting
data <- data %>%
          remove_rownames() %>%
          column_to_rownames('Mouse.Number') %>% 
          select(-Sex, -DOB, -Tissue.Collection.Date, -`Batch.GDF11/MSTN`) %>%
          select(-grep('LV|Fraction|Stroke|Cardiac.Output', colnames(.))) %>%
          rename(body.weight    = Body.Weight.g,
                 heart.weight   = Heart.Weight.mg,
                 tibia.length   = Tibia.Length.mm,
                 hw.bw.percent  = `HW/BW.%`,
                 hw.bw.ratio    = `HW/BW.mg/g`,
                 hw.tl.ratio    = `HW/TL.mg/mm`,
                 GDF11          = `GDF11.ng/mL`,
                 MSTN           = `MSTN.ng/mL`,
                 wall.thickness = `Wall.Thickness.µm`,
                 cardiomyocyte.cross.sectional.area = `Cardiomyocyte.Cross.Sectional.Area.µm2`,
                 percent.fibrosis = `Percent.Fibrosis.% area`)








### Normalize data
norm <- log(data)


### Temporary data frame for regression
temp <- cbind(model.matrix(~sex, data = annot.sample)[,-1,drop = FALSE], norm)

### Regress out body weight in heart weight data
fit  <- lm(heart.weight ~ sexM + body.weight, data = temp, na.action = na.exclude)
coef <- coefficients(fit)
res  <- residuals(fit)
hw.adj.bw <- t(coef[1] + (coef[2] %*% temp[,names(coef)[2]]) + res)

### Regress out tibia length in heart weight data
fit  <- lm(heart.weight ~ sexM + tibia.length, data = temp, na.action = na.exclude)
coef <- coefficients(fit)
res  <- residuals(fit)
hw.adj.tl <- t(coef[1] + (coef[2] %*% temp[,names(coef)[2]]) + res)



### Combined adjusted values to normalized data
norm <- cbind(norm, hw.adj.bw, hw.adj.tl)






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





### Covar info
covar.info <- data.frame(sample.column = c('sex', 'batch'),
                         covar.column  = c('sex', 'batch'),
                         display.name  = c('Sex', 'Batch'),
                         interactive   = c(TRUE, FALSE),
                         primary       = c('sex', NA),
                         lod.peaks     = c('sex_int', NA))
















### Creating annot.phenotype dataframe
annot.phenotype <- data.frame(data.name      = c(colnames(annot.sample), colnames(norm)),
                              short.name     = c(colnames(annot.sample), colnames(norm)),
                              R.name         = c(colnames(annot.sample), colnames(norm)),
                              description    = c('Mouse identifier', 'Sex of mouse: Female (F) or Male (M)', 'Mass spectromery batch group for GDF11 and MSTN', 'Mouse data-of-birth', 'Tissue collection date',
                                                 'Body weight measured in grams', 'Heart weight measured in milligrams', ' Tibia length measured in millimeters', 
                                                 '(Heart weight / 100) / (Body weight / 100)', 'Heart weight to body weight ratio',
                                                 'Heart weight to tibia length ratio', 'GDF11 levels', 'MSTN levels', 'Thickness of heart wall measured in micrometer',
                                                 'Cardiomyocyte cross sectional area', 'Fibrosis percentage', 'Heart weight conditioned on body weight',
                                                 'Heart weight conditioned on tibia length'),
                              units          = c(rep(NA, ncol(annot.sample)), 'g', 'mg', 'mm', '%', 'mg/g', 'mg/mm', 'ng/mL', 'ng/mL', 'µm', 'µm2', '% area', NA, NA),
                              category       = c(rep('Demographic', ncol(annot.sample)), rep('Phenotype', ncol(norm))),
                              R.category     = c(rep('Demographic', ncol(annot.sample)), rep('Phenotype', ncol(norm))),
                              is.id          = c(TRUE, rep(FALSE, ncol(annot.sample) + ncol(norm) - 1)),
                              is.numeric     = c(rep(FALSE, ncol(annot.sample)), rep(TRUE, ncol(norm))),
                              is.date        = c(FALSE, FALSE , FALSE, TRUE, TRUE, rep(FALSE, ncol(norm))),
                              is.factor      = c(FALSE, TRUE, TRUE, FALSE , FALSE, rep(TRUE, ncol(norm))),
                              factor.levels  = c(NA, 'F:M', '1:2', NA, NA, rep(NA, ncol(norm))),
                              is.covar       = c(FALSE, TRUE, TRUE, FALSE , FALSE, rep(FALSE, ncol(norm))),
                              is.pheno       = c(rep(FALSE, ncol(annot.sample)), rep(TRUE, ncol(norm))),
                              is.derived     = c(rep(FALSE, ncol(annot.sample)), FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, rep(FALSE, ncol(norm) - 8), TRUE, TRUE),
                              omit           = FALSE,
                              use.covar      = c(rep(NA, ncol(annot.sample)), rep('sex', 6), 'sex:batch', 'sex:batch', rep('sex', 5)),
                              transformation = c(rep(NA, ncol(annot.sample)), rep('log', ncol(norm))),
                              adj.pheno      = c(rep(NA, ncol(annot.sample)), rep(NA, ncol(norm) - 2), 'body.weight', 'tibia.length'))

















### Format for QTL viewer
dataset.heart.phenotype <- list(annot.phenotype = as_tibble(annot.phenotype),
                                annot.samples   = as_tibble(annot.sample),
                                covar.matrix    = covar.matrix,
                                covar.info      = as_tibble(covar.info),
                                data            = list(raw  = as.matrix(data),
                                                       norm = as.matrix(norm),
                                                       rz   = as.matrix(rz)),
                                datatype        = 'phenotype',
                                display.name    = 'Heart physiological phenotypes',
                                lod.peaks       = list())







### Save
rm(list = ls()[!grepl('dataset[.]|genoprobs|K|map|markers', ls())])
save.image(file = 'pazdro_heart_phenotype_viewer_v1.Rdata')                              
                              
