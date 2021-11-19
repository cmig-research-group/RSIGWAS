#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

infile = args[1]

require('dplyr')
require('ggplot2')
require('reshape2')
require('lme4')
require('R.matlab')
require('lmerTest')

# Assemble data ---- release 2.0.1
abcd <- readRDS('~/Phenotypes/nda3.0.Rds')
abcd$IID <- abcd$subjectid
pcs <- read.table('~/ABCD/PCs/plink2.eigenvec', header=T)
ethnicity <- read.table('~/ABCD/PCs/race_ethnicity.tsv', header=T)

# Get joint 
abcd <- inner_join(abcd, pcs)
abcd <- inner_join(abcd, ethnicity)

covar_names <- c('age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','mri_info_device.serial.number','smri_vol_subcort.aseg_intracranialvolume','genetic_ancestry_factor_european')
abcd.sub <- abcd[c(covar_names, 'eventname','rel_family_id')]
abcd.sub$IID <- abcd$subjectid

mat = readMat(infile)

df <- data.frame()
for (idx in 1:dim(mat$vol)[2]) {
  cat(paste0(idx, '\r\n'))
  test <- data.frame(IID = unlist(mat$subjid), event = unlist(mat$eventvec), y = mat$vol[,idx])
  test$eventname <- 'baseline_year_1_arm_1'
  test$eventname[test$event == '2year'] <- '2_year_follow_up_y_arm_1'
  test.geno <- data.frame(IID = unlist(mat$genoid), x = mat$geno[,idx])
  tmp <- left_join(test, abcd.sub)
  tmp <- left_join(tmp, test.geno)
  tmp$EUR <- tmp$genetic_ancestry_factor_european > 0.8
  tmp <- subset(tmp, !is.na(tmp$EUR))
  tmp <- subset(tmp, !duplicated(tmp$IID))
  h1 <- as.formula(paste0('y ~ ', paste0(covar_names[1:(length(covar_names)-1)], collapse='+'), ' + x + (1|rel_family_id) + ( 1 | EUR)'))
  fit <- lmer(h1, data=tmp)
  df_tmp <- summary(fit)$coefficients[dim(summary(fit)$coefficients)[1],]
  df <- rbind(df, df_tmp)
}
names(df) <- names(df_tmp)

writeMat(paste0(infile, '.R.results.mat'), df=df)





