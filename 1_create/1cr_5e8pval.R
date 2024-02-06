################################################################################
#	R-script:		  1cr_5e8pval.R
#	Project:	  	PROECMR : Protein and EC Risk - Mendelian randomisation
#
#	Data used:		1_initialdata/proteins
#
#	Data created:	2_derivedata/proteins_p5e8
#
# Purpose:  		cr dataset with variants p<5e-8 for each protein
################################################################################

## set file path
rm(list = ls(all.names = TRUE))
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/1_initialdata/proteins"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_p5e8"

## cr list of file names
setwd(INPUT)
file_list <- list.files()

## loop through each file
for (file_name in file_list) {
  # load data
  setwd(PROECMRinitial)
  data_olk=read.table(file_name, sep = "", header = TRUE)

  # convert logpval to numeric
  data_olk$logpval <-  c(as.numeric(data_olk$LOG10P))
  
  # subset variants >-log10(5e-8)
  data_olk_sigp <- subset(data_olk, logpval > -log10(5e-8))
  
  # cr csv file
  setwd(OUTPUT)
  write.table(data_olk_sigp,  file = paste0(file_name, "_p5e8.txt"))
}


