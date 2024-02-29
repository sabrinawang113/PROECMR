
rm(list=ls())

# environment ====
library(dplyr)
library(tidyr)
library(data.table)
library(TwoSampleMR)
library(ieugwasr)

## file paths ====
OUTPUT = "C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/2_derivedata/mr_il6/"

## data_exposure ====
data_exposure_all <- read.table("C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/1_initialdata/UKBPPP_Sun/IL6_P05231_OID20563", sep = "", header = TRUE)
data_exposure_all$LOG10P <-  as.numeric(data_exposure_all$LOG10P)

data_exposure <- subset(data_exposure_all, LOG10P > -log10(5e-8))
data_exposure$BETA <-  as.numeric(data_exposure$BETA)
data_exposure$SE <-  as.numeric(data_exposure$SE)
data_exposure$A1FREQ <-  as.numeric(data_exposure$A1FREQ)
label <- data_exposure$phenotype[[1]]

data_exposure <- format_data(data_exposure, type="exposure", snps = NULL,  header = TRUE, 
                           phenotype_col = "phenotype", 
                           id_col = "phenotype", 
                           snp_col = "rsid", 
                           beta_col = "BETA", 
                           se_col = "SE", 
                           pval_col = "LOG10P", log_pval = TRUE,
                           eaf_col = "A1FREQ", 
                           effect_allele_col = "ALLELE1", 
                           other_allele_col = "ALLELE0", 
                           chr_col = "CHROM", 
                           pos_col = "GENPOS",
                           samplesize_col = "N")
data_exposure$id.exposure <- label
write.table(data_exposure, paste0(OUTPUT, "data_exposure.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## data_clump ====
colnames(data_exposure)[colnames(data_exposure) == "pval"] <- "oldpval"
colnames(data_exposure)[colnames(data_exposure) == "SNP"] <- "rsid"
colnames(data_exposure)[colnames(data_exposure) == "pval.exposure"] <- "pval"

data_clumped <- ld_clump(dat = data_exposure,
                         clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
                         pop = "EUR",
                         plink_bin = genetics.binaRies::get_plink_binary())
colnames(data_clumped)[colnames(data_clumped) == "rsid"] <- "SNP"
colnames(data_clumped)[colnames(data_clumped) == "pval"] <- "pval.exposure"

data_clumped$f_stats <- (data_clumped$b / data_clumped$se)^2 
write.table(data_clumped, paste0(OUTPUT, "data_clump_il6.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## data_outcome ====
data_outcome <- read.table("C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/1_initialdata/ECAC/ECAC2018_noUKBB_strata.txt", sep = "", header = TRUE)
data_outcome$phenotype <- "ECAC_noUKBB"
data_outcome <- format_data(data_outcome, type="outcome", snps = data_clumped$SNP,  header = TRUE, 
                                  phenotype_col = "phenotype", 
                                  snp_col = "SNPID", 
                                  beta_col = "BETA", 
                                  se_col = "SE", 
                                  eaf_col = "MEAN_EAF", 
                                  effect_allele_col = "EA", 
                                  other_allele_col = "OA", 
                                  pval_col = "PVALUE", 
                                  chr_col = "CHR", 
                                  pos_col = "POS")
data_outcome$id.outcome <- "ECAC_noUKBB"
write.table(data_outcome, paste0(OUTPUT, "data_outcome.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## harmonise ====
data_harmonise <- harmonise_data(data_clumped, data_outcome, action=2)
write.table(data_harmonise, paste0(OUTPUT, "data_harmonise.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# steiger filtering ====
data_harmonise$ncase.outcome <- 12270
data_harmonise$ncontrol.outcome <- 46126
data_harmonise$prevalence.outcome <- 0.000158
data_harmonise$units.outcome <- "log odds"
data_steiger <- steiger_filtering(data_harmonise)
write.table(data_steiger, paste0(OUTPUT, "data_steiger.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## mr ====
methods <- mr_method_list()
methods_heterogeneity <- subset(methods, heterogeneity_test == TRUE)$obj
methods <- c(
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_ivw_fe",
  "mr_weighted_mode")
data_mr <- mr(data_harmonise, method_list = methods)
write.table(data_mr, paste0(OUTPUT, "data_mr.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## sensitivity anlayses ====
data_sensitivity <- mr_singlesnp(data_harmonise)
write.table(data_sensitivity, paste0(OUTPUT, "data_mr-singlesnp.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
plot_sensitivity <- mr_forest_plot(data_sensitivity)
pdf(paste0(OUTPUT, "plot_mr-singlesnp.pdf"))
for (i in 1:length(plot_sensitivity)) {
  print(plot_sensitivity[[i]])
}
dev.off()

data_sensitivity <- mr_leaveoneout(data_harmonise, method = mr_ivw_fe)
write.table(data_sensitivity, paste0(OUTPUT, "data_mr-leaveoneout.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
plot_sensitivity <- mr_leaveoneout_plot(data_sensitivity)
pdf(paste0(OUTPUT, "plot_mr-leaveoneout.pdf"))
for (i in 1:length(plot_sensitivity)) {
  print(plot_sensitivity[[i]])
}
dev.off()

data_harmonise <- fread(paste0(OUTPUT, "data_harmonise.txt"))
data_mr <- fread(paste0(OUTPUT, "data_mr.txt"))
plot_sensitivity <- mr_scatter_plot(data_mr, data_harmonise)
pdf(paste0(OUTPUT, "plot_mr-scatter.pdf"))
for (i in 1:length(plot_sensitivity)) {
  print(plot_sensitivity[[i]])
}
dev.off()