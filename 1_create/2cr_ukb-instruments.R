################################################################################
# R-script: 2cr_ukb-instruments.R
# Project: PROECMR : Protein and EC Risk - Mendelian randomisation
#
# Purpose: select cis (or trans) variants as iv for each protein & LD clump
#
################################################################################

## environment
rm(list = ls(all.names = TRUE))
#library(devtools)
#install_github("MRCIEU/ieugwasr")
library(ieugwasr)

## set file path
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_5e8"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_instrument_1mb"

## data: file gene position to determine cis variants
data_genepos=read.table("~/_PROJECTS/Proj_PROECMR_osiris/1_initialdata/olink_protein_map_1.5k_v1.txt" , sep = "\t", header = TRUE)

## selects cis snps then clump
# cr list of file names
setwd(INPUT)
file_list <- list.files()

file_list <- file_list[c(1:4,9:12)]

# loop through each file
for (file_name in file_list) {
  # load data
  setwd(INPUT)
  data_olk=read.table(file_name,  header = TRUE)
  
  # obtain chrom & gene start and end 
  olk_name <- substr(data_olk$phenotype, 1, regexpr("_", data_olk$phenotype) - 1)
  olk_name <- olk_name[1]
  gene_chr <- subset(data_genepos$chr, data_genepos$Assay == olk_name)
  gene_start <- subset(data_genepos$gene_start, data_genepos$Assay == olk_name) - 1000000
  gene_end <- subset(data_genepos$gene_end, data_genepos$Assay == olk_name) + 1000000
  
  ## subset cis variants
  data_olk$GENPOS <-  c(as.numeric(data_olk$GENPOS))
  data_olk_cis <- subset(data_olk, CHROM==gene_chr & GENPOS>=gene_start & GENPOS<=gene_end)
  
  ## clumping
  # make sure col names = rsid & pval
  # if no cis variant, clump trans variants
  if (nrow(data_olk_cis) == 0) {
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
  # if yes cis variants, clump cis variants
  if (nrow(data_olk_cis) != 0) {
    data_olk_cis$pval <-  c(10^(-data_olk_cis$logpval)) # -logpval -> pval
    data_clump <- ld_clump(dat = data_olk_cis,
                           clump_kb = 10000, clump_r2 = 0.001,
                           pop = "EUR",
                           access_token = NULL,
                           bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
                           plink_bin = genetics.binaRies::get_plink_binary())
    # export data
    setwd(OUTPUT)
    write.table(data_clump,  file = paste0(file_name, "_cis.txt"))
  }
}

################################################################################

## select IL6R missense variant rs2228145 
file_name <- c("IL6R_P08887_OID20385_p5e8.txt")
# load data
setwd(INPUT)
data_olk=read.table(file_name,  header = TRUE)
data_olk$pval <-  c(10^(-data_olk$logpval)) # -logpval -> pval
data_olk$id <-  "IL6R" # -logpval -> pval
data_clump <- subset(data_olk, data_olk$rsid == "rs2228145")
# export data
setwd(OUTPUT)
write.table(data_clump,  file = "IL6R_rs2228145_cis.txt")

################################################################################

# list OLK files
setwd(OUTPUT)
file_list <- list.files()
ivsnps_clump <- data.frame()


#file_name <- file_list[5]
#data_olk=read.table(file_name,  header = TRUE)
#data_olk <- subset(data_olk, data_olk$logpval>8)
#write.table(data_olk,  file = file_name)


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
write.table(ivsnps_clump,  file = "list_proteins_instrument.txt")



setwd("/home/wangs/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_instrument_1mb")
file_list <- list.files()
ivsnps_clump <- data.frame()
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
ivsnps_clump %>%
  group_by(phenotype) %>%
  summarise(count = n())
