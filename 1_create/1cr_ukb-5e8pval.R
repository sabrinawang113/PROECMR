################################################################################
# R-script: 1cr_ukb-5e8pval.R
# Project: PROECMR : Protein and EC Risk - Mendelian randomisation
#
# Data used: 1_initialdata/proteins
#
# Data created: 2_derivedata/proteins_p5e8
#
# Purpose: cr dataset with variants p<5e-8 for each protein
#
################################################################################

## set file path
rm(list = ls(all.names = TRUE))
INPUT <- "/home/wangs/_PROJECTS/Proj_PROECMR_osiris/1_initialdata/proteins_combined"
OUTPUT <- "/home/wangs/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_5e8"


## cr list of file names
setwd(INPUT)
file_list <- list.files()
file_list <- file_list[3:12]


## loop through each file
for (file_name in file_list) {
  # load data
  print(paste0("loading:", file_name))
  data_olk=read.table(file_name, header = TRUE)

  # convert logpval to numeric
  print(paste0("converting logp:", file_name))
  data_olk$logpval <-  c(as.numeric(data_olk$LOG10P))
  
  # subset variants >-log10(5e-8)
  print(paste0("subsetting 5e-8:", file_name))
  data_olk_sigp <- subset(data_olk, logpval > -log10(5e-8))
  
  # cr csv file
  print(paste0("saving:", file_name))
  write.table(data_olk_sigp,  file = paste0("/home/wangs/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_5e8/", file_name, "_p5e8.txt"))
}