rm(list = ls(all.names = TRUE))

## load packages
library(dplyr)
library(MendelianRandomization)

## set file path
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_chargeukb"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/5_outputs"

## load ECAC data
data_ecac = read.table("/home/wangs/_PROJECTS/Proj_PROECMR_osiris/1_initialdata/ECAC2018_noUKBB_strata.txt", sep = "", header = TRUE)
colnames(data_ecac)[colnames(data_ecac) == "POS"] <- "POS19"
colnames(data_ecac)[colnames(data_ecac) == "CHR"] <- "CHROM"

## empty dataframe
mr_inputs <- data.frame()
mr_results <- data.frame()

setwd(INPUT)
file_list <- list.files()
print(file_list)

## loop MR analysis
for (file_name in file_list) {
  # load data
  setwd(INPUT)
  data_chargeukb=read.table(file_name,  header = TRUE)
  
  # macro protein name
  label <- sub("\\.txt$", "", file_name)
  
  ## select snps
  snps <- data_chargeukb$SNP
  snp_ecac <- subset(data_ecac, SNPID %in% c(snps))
  snp_chargeukb <- data_chargeukb
  
  # rename ECAC rsid, effectallele, refallele, eafreq, beta, se
  colnames(snp_ecac)[colnames(snp_ecac) == "SNPID"] <- "rsid"
  colnames(snp_ecac)[colnames(snp_ecac) == "EA"] <- "ecac_effectallele"
  colnames(snp_ecac)[colnames(snp_ecac) == "OA"] <- "ecac_refallele"
  colnames(snp_ecac)[colnames(snp_ecac) == "MEAN_EAF"] <- "ecac_eafreq"
  colnames(snp_ecac)[colnames(snp_ecac) == "BETA"] <- "ecac_beta"
  colnames(snp_ecac)[colnames(snp_ecac) == "SE"] <- "ecac_se"
  colnames(snp_ecac)[colnames(snp_ecac) == "PVALUE"] <- "ecac_pval"
  snp_ecac$olk <- label 
  # rename CHARGE rsid, effectallele, refallele, eafreq, beta, se
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "SNP"] <- "rsid"
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "effect_allele.exposure"] <- "olk_effectallele"
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "other_allele.exposure"] <- "olk_refallele"
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "beta.exposure"] <- "olk_beta"
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "se.exposure"] <- "olk_se"
  colnames(snp_chargeukb)[colnames(snp_chargeukb) == "pval.exposure"] <- "olk_pval"
  snp_chargeukb$olk <- label
  
  ## align effect alleles and change outcome beta sign if needed
  snp_merge <- data.frame(merge(snp_chargeukb, snp_ecac, by = "rsid"))
  snp_merge <- snp_merge %>%
    mutate(ecac_beta_correct = ifelse(olk_effectallele == ecac_effectallele, ecac_beta, ecac_beta * -1))
  
  ## invert beta so that estimates represent decreased crp
  snp_merge$olk_beta <- (snp_merge$olk_beta)*-1
  
  ## set up MR object
  bx.all = snp_merge$olk_beta
  bxse.all = snp_merge$olk_se
  by.all = snp_merge$ecac_beta_correct
  byse.all = snp_merge$ecac_se
  MRObject=mr_input(bx=bx.all,bxse = bxse.all, by=by.all, byse=byse.all)
  MRObject
  
  # post results to df
  mr_object <- data.frame(olk=label, rsid=snp_merge$rsid, bx=bx.all,bxse = bxse.all, by=by.all, byse=byse.all)
  mr_inputs <- rbind(mr_inputs, mr_object)
  
  
  ## run MR
  mr_ivw(MRObject,model = "fixed")
  mr_estimates <-mr_ivw(MRObject,model = "fixed")
  
  ivw=c(mr_estimates@Model,
        label,
        mr_estimates@Estimate,
        mr_estimates@StdError,
        mr_estimates@CILower,
        mr_estimates@CIUpper,
        mr_estimates@Pvalue,
        mr_estimates@Heter.Stat,
        mr_estimates@Fstat)
  ivw <- data.frame(
    Model = ivw[[1]],
    Exposure = ivw[[2]],
    Estimate = ivw[[3]],
    StdError = ivw[[4]],
    CILower = ivw[[5]],
    CIUpper = ivw[[6]],
    Pvalue = ivw[[7]],
    HetStatQ = ivw[[8]],
    HetStatPval = ivw[[9]],
    Fstat = ivw[[10]]
  )
  
  # post results to df
  mr_results <- rbind(mr_results, ivw)
}

## export results
setwd(OUTPUT)
write.table(mr_inputs,  file = "mr_inputs_charge.txt")
write.table(mr_results,  file = "mr_results_charge.txt")
