################################################################################
#	R-script:		  4cr_chargeukb.R
#	Project:	  	PROECMR : Protein and EC Risk - Mendelian randomisation
#
#	Data used:		1_initialdata/CHARGEUKB_GCST90029070_buildGRCh37.tsv.gz
#               
#	Data created:	2_derivedata/proteins_chargeukb
#
# Purpose:  		select cis variants as iv for IL6R and CRP from CHARGEUKB GWAS
#               LD clump to select independent snps
#               r<0.01 for cis variants
################################################################################

## set file path
rm(list = ls(all.names = TRUE))
INPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/1_initialdata"
OUTPUT <- "~/_PROJECTS/Proj_PROECMR_osiris/2_derivedata/proteins_chargeukb"

library(tidyverse)
library(TwoSampleMR)
library(dplyr)
library(MendelianRandomization)

## data: CHARGEUKB GWAS
setwd(INPUT)
crp  <- data.table::fread("CHARGEUKB_GCST90029070_buildGRCh37.tsv.gz") 

# gene region for crp and il6r
crp_start <- 159682079 
crp_end <- 159684379 

IL6R_start <- 154377819 
IL6R_end <- 154441926 

# ld clump
crp_clumped <- crp %>%
  filter(chromosome == 1 & (base_pair_location > crp_start) & (base_pair_location < crp_end)) %>%
  dplyr::select(SNP = variant_id, effect_allele.exposure = effect_allele, other_allele.exposure = other_allele, beta.exposure = beta, se.exposure = standard_error, pval.exposure = p_value) %>%
  as_tibble() %>%
  mutate(eaf.exposure = NA) %>%
  mutate(exposure ="CRP", id.exposure = "CRP") %>%
  mutate(pval.exposure = as.numeric(pval.exposure)) %>%
  filter(pval.exposure < 5e-8)  %>%
  clump_data(clump_r2 = 0.01)


il6r_clumped <- crp %>%
  filter(chromosome == 1 & (base_pair_location > IL6R_start) & (base_pair_location < IL6R_end)) %>%
  select(SNP = variant_id, effect_allele.exposure = effect_allele, other_allele.exposure = other_allele, beta.exposure = beta, se.exposure = standard_error, pval.exposure = p_value) %>%
  as_tibble() %>%
  mutate(eaf.exposure = NA) %>%
  mutate(exposure ="IL-6", id.exposure = "IL-6") %>%
  mutate(pval.exposure = as.numeric(pval.exposure)) %>%
  filter(pval.exposure < 5e-8)  %>%
  clump_data(clump_r2 = 0.01)


# Save results
setwd(OUTPUT)
write.table(crp_clumped,  file = "CHARGEUKB_crp_cis.txt")
write.table(il6r_clumped,  file = "CHARGEUKB_il6r_cis.txt")

