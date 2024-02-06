################################################################################
#	R-script:		  3an_mr_ivw_sensitivityanalyses.R
#	Project:	  	PROECMR : Protein and EC Risk - Mendelian randomisation
#
#	Data used:		Proteins: 2_derivedata/proteins_instrument_trans
#               Cancer:   ECAC2018_noUKBB_strata.txt
#
#	Data created:	sensitivityanalyses/mr_rsid_olk.txt
#               sensitivityanalyses/mr_rsid_ecac.txt
#               sensitivityanalyses/mr_inputs.txt
#               sensitivityanalyses/mr_results.txt
#               sensitivityanalyses/mr_plot/
#
# Purpose:  		sensitivity analysis perform two-sample MR with cis/trans snps
#               
################################################################################

## load packages
library(dplyr)
library(MendelianRandomization)

## set file path
rm(list = ls(all.names = TRUE))
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_instrument_trans"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/5_outputs/sensitivityanalyses"

## load ECAC data
data_ecac = read.table("/home/wangs/_PROJECTS/Proj_PROECMR_osiris/1_initialdata/ECAC2018_noUKBB_strata.txt", sep = "", header = TRUE)
colnames(data_ecac)[colnames(data_ecac) == "POS"] <- "POS19"
colnames(data_ecac)[colnames(data_ecac) == "CHR"] <- "CHROM"

################################################################################

## empty dataframe
mr_rsid_olk <- data.frame()
mr_rsid_ecac <- data.frame()
mr_inputs <- data.frame()
mr_results <- data.frame()

setwd(INPUT)
file_list <- list.files()
print(file_list)

## loop MR analysis
for (file_name in file_list) {
  # load data
  setwd(INPUT)
  data_olk=read.table(file_name,  header = TRUE)
  
  # macro olk protein name
  olk_name <- substr(file_name, 1, regexpr("_", file_name) - 1)
  
  ## select snps
  chrom <- data_olk$CHROM
  pos19 <- data_olk$POS19
  snp_ecac <- subset(data_ecac, CHROM %in% c(chrom))
  snp_ecac <- subset(snp_ecac, POS19 %in% c(pos19))
  snp_olk <- data_olk
  
  # rename ECAC rsid, effectallele, refallele, eafreq, beta, se
  colnames(snp_ecac)[colnames(snp_ecac) == "SNPID"] <- "rsid"
  colnames(snp_ecac)[colnames(snp_ecac) == "EA"] <- "ecac_effectallele"
  colnames(snp_ecac)[colnames(snp_ecac) == "OA"] <- "ecac_refallele"
  colnames(snp_ecac)[colnames(snp_ecac) == "MEAN_EAF"] <- "ecac_eafreq"
  colnames(snp_ecac)[colnames(snp_ecac) == "BETA"] <- "ecac_beta"
  colnames(snp_ecac)[colnames(snp_ecac) == "SE"] <- "ecac_se"
  colnames(snp_ecac)[colnames(snp_ecac) == "PVALUE"] <- "ecac_pval"
  snp_ecac$olk <- olk_name 
  # rename OLK rsid, effectallele, refallele, eafreq, beta, se
  colnames(snp_olk)[colnames(snp_olk) == "ALLELE1"] <- "olk_effectallele"
  colnames(snp_olk)[colnames(snp_olk) == "ALLELE0"] <- "olk_refallele"
  colnames(snp_olk)[colnames(snp_olk) == "A1FREQ"] <- "olk_eafreq"
  colnames(snp_olk)[colnames(snp_olk) == "BETA"] <- "olk_beta"
  colnames(snp_olk)[colnames(snp_olk) == "SE"] <- "olk_se"
  colnames(snp_olk)[colnames(snp_olk) == "pval"] <- "olk_pval"
  snp_olk$olk <- olk_name
  
  # post results to df
  mr_rsid_olk <- rbind(mr_rsid_olk, snp_olk)
  mr_rsid_ecac <- rbind(mr_rsid_ecac, snp_ecac)
  
  ## align effect alleles and change outcome beta sign if needed
  snp_merge <- data.frame(merge(snp_olk, snp_ecac, by = "rsid"))
  snp_merge <- snp_merge %>%
    mutate(ecac_beta_correct = ifelse(olk_effectallele == ecac_effectallele, ecac_beta, ecac_beta * -1))
  snp_merge <- snp_merge %>%
    mutate(ecac_eafreq_correct = ifelse(olk_effectallele == ecac_effectallele, ecac_eafreq, 1-ecac_eafreq))
  
  ## set up MR object
  bx.all = snp_merge$olk_beta
  bxse.all = snp_merge$olk_se
  by.all = snp_merge$ecac_beta_correct
  byse.all = snp_merge$ecac_se
  MRObject=mr_input(bx=bx.all,bxse = bxse.all, by=by.all, byse=byse.all)
  MRObject
  
  # post results to df
  mr_object <- data.frame(olk=olk_name, rsid=snp_merge$rsid, bx=bx.all,bxse = bxse.all, by=by.all, byse=byse.all)
  mr_inputs <- rbind(mr_inputs, mr_object)
  
  
## run MR: ivw
  mr_ivw(MRObject,model = "random")
  mr_estimates <-mr_ivw(MRObject,model = "random")
  
  result=c(mr_estimates@Model,
        olk_name,
        mr_estimates@Estimate,
        mr_estimates@StdError,
        mr_estimates@CILower,
        mr_estimates@CIUpper,
        mr_estimates@Pvalue)
  ivw <- data.frame(
    Method = "IVW",
    Exposure = result[[2]],
    Estimate = result[[3]],
    StdError = result[[4]],
    CILower = result[[5]],
    CIUpper = result[[6]],
    Pvalue = result[[7]]
  )

## run MR: egger
mr_estimates <-mr_egger(MRObject)

result=c(mr_estimates@Model,
      olk_name,
      mr_estimates@Estimate,
      mr_estimates@StdError.Est,
      mr_estimates@CILower.Est,
      mr_estimates@CIUpper.Est,
      mr_estimates@Pvalue.Est)
egger <- data.frame(
  Method = "Egger",
  Exposure = result[[2]],
  Estimate = result[[3]],
  StdError = result[[4]],
  CILower = result[[5]],
  CIUpper = result[[6]],
  Pvalue = result[[7]]
)

## run MR: mr_median
mr_estimates <-mr_median(MRObject)

result=c(mr_estimates@Type,
      olk_name,
      mr_estimates@Estimate,
      mr_estimates@StdError,
      mr_estimates@CILower,
      mr_estimates@CIUpper,
      mr_estimates@Pvalue)
weightedmed <- data.frame(
  Method = "WeightedMed",
  Exposure = result[[2]],
  Estimate = result[[3]],
  StdError = result[[4]],
  CILower = result[[5]],
  CIUpper = result[[6]],
  Pvalue = result[[7]]
)

# post results to df
mr_results <- rbind(mr_results, ivw, egger, weightedmed)

## plot
MRAllObject_main <- mr_allmethods(MRObject, method = "main")
pdfname <- paste(OUTPUT, "/mr_plot_",olk_name, ".pdf", sep = "")
pdf(pdfname)
mr_plot <- mr_plot(MRAllObject_main)
print(mr_plot)
dev.off()

}

## export results
setwd(OUTPUT)
write.table(mr_rsid_olk,  file = "mr_rsid_olk.txt")
write.table(mr_rsid_ecac,  file = "mr_rsid_ecac.txt")
write.table(mr_inputs,  file = "mr_inputs.txt")
write.table(mr_results,  file = "mr_results.txt")
