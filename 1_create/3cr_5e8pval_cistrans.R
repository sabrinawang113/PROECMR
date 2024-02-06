################################################################################
#	R-script:		  3cr_5e8pval_cistrans.R
#	Project:	  	PROECMR : Protein and EC Risk - Mendelian randomisation
#
#	Data used:		2_derivedata/proteins_p5e8
#
#	Data created:	2_derivedata/proteins_instrument_trans
#
# Purpose:  		sensitivity analysis:
#               select cis/trans variants as iv for each protein
#               clump to select independent snps
################################################################################

## set file path
rm(list = ls(all.names = TRUE))
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_p5e8"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_instrument_trans"

## install and load packages
#library(devtools)
#install_github("MRCIEU/ieugwasr")
library(ieugwasr)
library(TwoSampleMR)

################################################################################

# cr list of file names
setwd(INPUT)
file_list <- list.files()

# loop through each file
for (file_name in file_list) {
  # load data
  setwd(INPUT)
  data_olk=read.table(file_name,  header = TRUE)
  
  ## clumping
  # make sure col names = rsid & pval
  # clump trans variants
    data_olk$pval <-  c(10^(-data_olk$logpval)) # -logpval -> pval
    data_clump <- ld_clump(dat = data_olk,
                           clump_kb = 10000, clump_r2 = 0.001,
                           pop = "EUR",
                           access_token = NULL,
                           bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                           plink_bin = genetics.binaRies::get_plink_binary())
    # export data
    setwd(OUTPUT)
    write.table(data_clump,  file = paste0(file_name, "_trans.txt"))
}


## run separately for IL6R (to overcome very small p=0 in R) ====
file_name <- file_list[[8]]
  # load data
  setwd(INPUT)
  data_olk=read.table(file_name,  header = TRUE)
  
  ## clumping
  # make sure col names = rsid & pval
  # clump trans variants
  data_olk$pval_actual <-  c(10^(-data_olk$logpval)) # -logpval -> actual pval
  data_olk$pval <-  c(10^((-data_olk$logpval)/20)) # -logpval -> pval for clumping (to overcome very small p=0 in R)
  data_clump <- ld_clump(dat = data_olk,
                         clump_kb = 10000, clump_r2 = 0.001,
                         pop = "EUR",
                         access_token = NULL,
                         bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                         plink_bin = genetics.binaRies::get_plink_binary())
  # export data
  setwd(OUTPUT)
  write.table(data_clump,  file = paste0(file_name, "_trans.txt"))

# format IL6R df so it fits with the rest of proteins  
  setwd("~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_instrument_trans")
  file_list <- list.files()
  file_name <- file_list[grepl("IL6R", file_list)]
  data_olk=read.table(file_name,  header = TRUE)
  data_olk <- subset(data_olk, select = -pval)
  colnames(data_olk)[colnames(data_olk) == "pval_actual"] <- "pval"  
  write.table(data_olk,  file =file_name)
  
################################################################################

## cr a df of snps that will be used as IV
# empty df
rm(list = ls(all.names = TRUE))
ivsnps_clump <- data.frame()

# list OLK files
setwd(OUTPUT)
file_list <- list.files()

# loop through OLK files
for (file_name in file_list) {
  # load data
  data_olk=read.table(file_name,  header = TRUE)
  
  # flag whether snps are cis or trans
  cis_trans <- sub(".*_(.*)", "\\1", file_name)
  
  ## cr a df of iv snps
  new_df <- data_olk[, c("phenotype", "rsid", "CHROM", "GENPOS", "POS19", "POS38", "ALLELE0", "ALLELE1",  "A1FREQ", "BETA", "SE", "logpval" )]
  new_df <- new_df %>%
    mutate(cis_trans = cis_trans)
  ivsnps_clump <- rbind(ivsnps_clump, new_df)
}

## export results
setwd("~/_PROJECTS/Proj_PROECMR_osiris/5_outputs")
write.table(ivsnps_clump,  file = "list_proteins_instrument_trans.txt")
