#正向
exposure_data<-fread("DrinksPerWeek.txt",header = T)
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

outcome_data<-fread("SmokingInitiation.txt",header = T)
merge<-merge(exp_dat,outcome_data,by.x = "SNP",by.y = "RSID")
write.csv(merge, file="merge.csv")
outcome_dat <-
  read_outcome_data(
    snps = exp_dat$SNP,
    filename = "merge.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "PVALUE"
  )

dat<-harmonise_data(exposure_dat = exp_dat,outcome_dat = outcome_dat)

dat <- dat[dat$mr_keep=="TRUE",]

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

#反向
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

outcome_data<-fread("DrinksPerWeek.txt",header = T)
merge<-merge(exp_dat,outcome_data,by.x = "SNP",by.y = "RSID")
write.csv(merge, file="merge.csv")
outcome_dat <-
  read_outcome_data(
    snps = exp_dat$SNP,
    filename = "merge.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    eaf_col = "AF",
    pval_col = "PVALUE"
  )

dat<-harmonise_data(exposure_dat = exp_dat,outcome_dat = outcome_dat)

dat <- dat[dat$mr_keep=="TRUE",]

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

