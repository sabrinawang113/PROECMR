rm(list=ls())

# environment ====
#devtools::install_github("jean997/cause@v1.2.0")
#devtools::install_github("explodecomputer/genetics.binaRies")
library("cause")
library(readr)
library(dplyr)
library(ieugwasr)
library(data.table)

# resource ====
# https://jean997.github.io/cause/ldl_cad.html

# data ====
data_exposure <- read.table("C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/1_initialdata/UKBPPP_Sun/IL6_P05231_OID20563", sep = "", header = TRUE)
data_exposure$BETA <-  as.numeric(data_exposure$BETA)
data_exposure$SE <-  as.numeric(data_exposure$SE)
data_outcome <- read.table("C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/1_initialdata/ECAC/ECAC2018_noUKBB_strata.txt", sep = "", header = TRUE)

# merge gwas for CAUSE
X <- cause::gwas_merge(data_exposure, data_outcome, snp_name_cols = c("rsid", "SNPID"), 
                beta_hat_cols = c("BETA", "BETA"), 
                se_cols = c("SE", "SE"), 
                A1_cols = c("ALLELE1", "EA"), 
                A2_cols = c("ALLELE0", "OA"))
X$pval1 = 2 * pnorm(-abs(X$beta_hat_1/X$seb1))
X$pval2 = 2 * pnorm(-abs(X$beta_hat_2/X$seb2))

# LD pruning ====
X_clump <- X %>%
  rename(rsid = snp,
         pval = pval1)
data_clump <- ieugwasr::ld_clump(dat = X_clump,
                     clump_r2 = 0.001,
                     clump_p = 5e-8,
                     plink_bin = genetics.binaRies::get_plink_binary(), 
                     pop = "EUR")
top_vars <- data_clump$rsid

# calculate nuisance parameters ====
set.seed(816)
varlist <- with(X, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(X, varlist)

params$rho
head(params$mix_grid)

# fit CAUSE ====
res <- cause(X=X, variants = top_vars, param_ests = params)
res$loos[[2]]
res$loos[[3]]
names(res)
res$elpd
summary(res, ci_size=0.95)

# save ====
save(X, params, data_clump, top_vars, file = "C:/Users/wangs/OneDrive - IARC/Projects_IARC/Proj_PROECMR/2_derivedata/mrcause_IL6_P05231_OID20563_combined.RData")

