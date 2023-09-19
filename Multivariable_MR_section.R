#This is only a example code
getwd()
setwd('E:/RWork')

library('plyr')
library('TwoSampleMR')
library('data.table')
library('tidyverse')
library('cause')
library('mrbayes')
library('MRPRESSO')
library('vcfR')
library('ieugwasr')
library('MendelianRandomization')
library('MVMR')###https://github.com/WSpiller/MVMR/blob/master/vignettes/MVMR.rmd
library('RColorBrewer')
library('cowplot')
library('gridExtra')
library('rjags') 
library('DescTools')
library('simex')
library('kableExtra')
library('RMediation')
library('tryx')
library('metaCCA')

#Drinks per week & Smoke MVMR
exposure_data_1 <- fread("DrinksPerWeek.txt", header = T)
exposure_data_P_1 <- subset(exposure_data_1, PVALUE < 5e-8)
write.csv(exposure_data_P_1, file = "exposure_data_Pexposure_1.csv")
C_Pexposure_1 <- read_exposure_data(
  filename = 'exposure_data_Pexposure_1.csv',
  clump = FALSE,
  sep= ",",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF",
  pval_col = "PVALUE"
)
exp_dat_1 <- clump_data(C_Pexposure_1, clump_r2 = 0.001, clump_kb = 10000)
X111 <- subset(C_Pexposure_1, SNP %in% exp_dat_1$SNP)

exposure_data_2<-fread("SmokingInitiation.txt",header = T)
exposure_data_P_2<-subset(exposure_data_2,PVALUE<5e-8)
write.csv(exposure_data_P_2, file="exposure_data_Pexposure_2.csv")
C_Pexposure_2 <- read_exposure_data(
  filename = 'exposure_data_Pexposure_2.csv',
  clump = FALSE,
  sep= ",",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF",
  pval_col = "PVALUE"
)
exp_dat_2 <-clump_data(C_Pexposure_2,clump_r2=0.001,clump_kb=10000)
X222<-subset(C_Pexposure_2,SNP %in% exp_dat_2$SNP)

outcome_data <- fread('finngen_R8_M13_DUPUTRYEN.txt',header = T)

SNP_list <- rbind.fill(X111, X222)

SNP_list_clump <- clump_data(SNP_list,clump_r2=0.001,clump_kb=10000)

SNP<-SNP_list_clump$SNP

SNP_u<-unique(SNP)

SNP<-as.data.frame(cbind(SNP_u,1))
colnames(SNP)<-c("SNP","nu")

Drink_X1<-merge(SNP,exposure_data_1,by.x = "SNP",by.y = "RSID")
Smoke_X2<-merge(SNP,exposure_data_2,by.x = "SNP",by.y = "RSID")

YG<-merge(SNP,outcome_data,by.x = "SNP",by.y = "rsids")

snps<-intersect(Drink_X1$SNP,Smoke_X2$SNP)
snps<-intersect(snps,YG$SNP)

d1<-format_data(
  exposure_data_1,
  type = "exposure",
  snps = snps,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF",
  pval_col = "PVALUE"
)

d1$exposure<-"Drink"
d1$id.exposure<-"Drink"
dh1 <- subset(d1, select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure))

d2<-format_data(
  exposure_data_2,
  type = "outcome",
  snps = snps,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "RSID",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "AF",
  pval_col = "PVALUE"
)

d2$outcome<-"Smoke"
d2$id.outcome<-"Smoke"

dh2 <- harmonise_data(d1, d2)

dh2 <- subset(dh2, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))

names(dh2) <- gsub("outcome", "exposure", names(dh2))

mv_exposures <- rbind(dh1, dh2)

d3<-format_data(
  outcome_data,
  type = "outcome",
  snps =  mv_exposures$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  eaf_col = "af_alt",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  pval_col = "pval",
)

d3$outcome<-"DD"
d3$id.outcome<-"DD"
dh3 <- subset(d3, select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome))

mv_DAT<-mv_harmonise_data(mv_exposures, dh3,  harmonise_strictness=2)
mv_DAT

res<-mv_multiple(mv_DAT)
res

exposure_beta<-as.data.frame(mv_DAT[["exposure_beta"]])
exposure_se<-as.data.frame(mv_DAT[["exposure_se"]])
exposure_se
outcome_beta<-mv_DAT[["outcome_beta"]]
outcome_se<-mv_DAT[["outcome_se"]]
SNP <- row.names(exposure_se)

mvmr<-as.data.frame(cbind(SNP,outcome_beta,outcome_se,exposure_beta$Drink,exposure_beta$Smoke,exposure_se$Drink,exposure_se$Smoke))
colnames(mvmr)<-c("SNP","betaYG","sebetaYG","betaX1","betaX2","sebetaX1","sebetaX2")
mvmr

mvmr$betaYG<-as.numeric(mvmr$betaYG)
mvmr$sebetaYG<-as.numeric(mvmr$sebetaYG)
mvmr$betaX1<-as.numeric(mvmr$betaX1)
mvmr$betaX2<-as.numeric(mvmr$betaX2)
mvmr$sebetaX1<-as.numeric(mvmr$sebetaX1)
mvmr$sebetaX2<-as.numeric(mvmr$sebetaX2)

mr_mvivw <- mr_mvivw(mr_mvinput(bx = cbind(mvmr$betaX1, mvmr$betaX2), bxse = cbind(mvmr$sebetaX1, mvmr$sebetaX2), by = mvmr$betaYG, mvmr$sebetaYG))
mr_mvegger <- mr_mvegger(mr_mvinput(bx = cbind(mvmr$betaX1, mvmr$betaX2), bxse = cbind(mvmr$sebetaX1, mvmr$sebetaX2), by = mvmr$betaYG, mvmr$sebetaYG))

mr_mvivw
mr_mvegger

head(mvmr)
bx <- as.matrix(mvmr[,c("betaX1", "betaX2")])
bxse <- as.matrix(mvmr[,c("sebetaX1", "sebetaX2")])
dat <- MendelianRandomization::mr_mvinput(bx = bx,
                                          bxse = bxse,
                                          by = mvmr$betaYG,
                                          byse = mvmr$sebetaYG,
                                          #snps = mvmr$SNP
)
dat <- mrmvinput_to_mvmr_format(dat)
head(dat)


strength_mvmr <- strength_mvmr(dat, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(dat, gencov=0)