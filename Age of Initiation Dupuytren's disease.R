###Age of Initiation Dupuytren's disease
getwd()
setwd('E:/RWork')

library('plyr')
library('TwoSampleMR')
library('data.table')
library('tidyverse')
library('MRPRESSO')
library('MendelianRandomization')
library('MVMR')
library('RColorBrewer')
library('cowplot')
library('gridExtra')
library('rjags') 
library('DescTools')
library('simex')
library('kableExtra')
library('tryx')
library('metaCCA')

exposure_data<-fread("AgeofInitiation.txt",header = T)
exposure_data_P<-subset(exposure_data,PVALUE<5e-8)
write.csv(exposure_data_P, file="exposure_data_Pexposure.csv")
C_Pexposure <-
  read_exposure_data(
    filename = "exposure_data_Pexposure.csv",
    sep = ",",
    snp_col = "RSID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "PVALUE",
    clump = FALSE
  )
exp_dat <-clump_data(C_Pexposure,clump_r2=0.001,clump_kb=10000)

outcome_data<-fread("finngen_R8_M13_DUPUTRYEN.txt",header = T)
merge<-merge(exp_dat,outcome_data,by.x = "SNP",by.y = "rsids")
write.csv(merge, file="merge.csv")
outcome_dat <-
  read_outcome_data(
    snps = exp_dat$SNP,
    filename = "merge.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    eaf_col = "af_alt"
  )

dat<-harmonise_data(exposure_dat = exp_dat,outcome_dat = outcome_dat)

dat <- dat[dat$mr_keep=="TRUE",]

dat$exposure <- "Age of Initiation"
dat$outcome <- "Palmar fascial fibromatosis [Dupuytren]s"

Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#Calculate Isq weighted and unweighted
I2<-c()

str(dat)

dat$samplesize.exposure <- 341427  
dat$samplesize.outcome <- 257738

#Rename required columns
dat $BetaXG<-dat $beta.exposure
dat $seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome
SNP<-dat $SNP
BXG = abs(BetaXG) 

# Calculate F statistics and I-squared statistics
F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted,SNP)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted","SNP")
write.csv(I2, file="regression_dilution_isq_weighted_afs.csv", row.names = FALSE)

#Evaluate primary MR results
res <-
  mr(
    dat,
    method_list = c(
      "mr_ivw",
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_raps",
      'mr_simple_median'
    )
  )

OR <-generate_odds_ratios(res)
OR

write.csv(OR, file="OR.csv")

#Heterogeneity
het <- mr_heterogeneity(dat)
write.csv( het, file="het.csv")

#Pleiotropy by MR-Egger 
pleio <- mr_pleiotropy_test(dat)
write.csv( pleio, file=" pleio.csv")

#Leave-one-out(Supplementary Figure 6. Leave-one-out analysis for age of initiation.)
single <- mr_leaveoneout(dat)
write.csv( single, file="single.csv")
mr_leaveoneout_plot( single)

#Scatter plot(Supplementary Figure 1. Scatter plot for age of initiation.)
mr_scatter_plot(res,dat)

#Funnel plot(Supplementary Figure 7. Funnel plot for age of initiation.)
mr_funnel_plot(res_single)

#MRPRESSO
PRESSO <-mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                   OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 5000,  
                   SignifThreshold = 0.05)

#Outliers output(If detected)
outliers<-dat$SNP[PRESSO[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]]
write.csv(outliers, file="A1outliers.csv")

#MRPRESSO-OR output
beta<-PRESSO[["Main MR results"]][["Causal Estimate"]][2]
se<-PRESSO[["Main MR results"]][["Sd"]][2]
OR = exp(beta)
UL_beta=beta + se * 1.96
LL_beta=beta - se * 1.96
LL_OR=exp(LL_beta)
UL_OR=exp(UL_beta)
Pvalue<-PRESSO[["Main MR results"]][["P-value"]][2]
presso_OR<-cbind(dat$exposure[1],dat$outcome[1],OR,LL_OR,UL_OR,Pvalue)
presso_OR

###Steiger filter(Example code)
dat <- steiger_filtering(dat)
dat <- dat %>% filter(steiger_dir == T)