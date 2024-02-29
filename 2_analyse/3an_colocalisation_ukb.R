rm(list=ls())

# environment ====
#if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
#remotes::install_github("MRCIEU/genetics.binaRies", force = F)
#remotes::install_github("explodecomputer/plinkbinr", force = F)
#remotes::install_github("chr1swallace/coloc@main", force = F)
#remotes::install_github("sjmgarnier/viridis", force = F)
library(genetics.binaRies)
library(plinkbinr)
library(coloc)
library(viridis)
library(data.table)
library(ieugwasr)
library(dplyr)
library(TwoSampleMR)
library(tidyverse)
source("/home/wangs/_PROJECTS/adiposity-proteins-endometrial/tools/my_coloc_chriswallace.R")


# file paths ====
DIRECTORY_OUT <- "/home/wangs/_PROJECTS/Proj_PROECMR_osiris/5_outputs/colocalization/"


# data_exposure ====
list_files_exposure <- list.files("/data/GWAS_data/work/UKB_PPP/cis-snps/european/500k", full.names = T)

CCL25 <- list_files_exposure[grepl("CCL25_", list_files_exposure)]
CLEC4G <- list_files_exposure[grepl("CLEC4G_", list_files_exposure)]
DNER <- list_files_exposure[grepl("DNER_", list_files_exposure)]
HGF <- list_files_exposure[grepl("HGF_", list_files_exposure)]
IL6R <- list_files_exposure[grepl("IL6R_", list_files_exposure)]
LILRB4 <- list_files_exposure[grepl("LILRB_", list_files_exposure)]
MASP1 <- list_files_exposure[grepl("MASP1_", list_files_exposure)]
PIK3AP1 <- list_files_exposure[grepl("PIK3AP1_", list_files_exposure)]
list_files_exposure <- c(CCL25,CLEC4G,DNER,HGF,IL6R,LILRB4,MASP1,PIK3AP1)

list_exposure <- lapply(list_files_exposure, fread, sep = " ", header = F)
label_exposure <- sub("/data/GWAS_data/work/UKB_PPP/cis-snps/european/500k/", "", list_files_exposure)
label_exposure <- sub(".gz", "", label_exposure)

## format ====
print("starting: for loop format exposure")
for (i in seq_along(list_exposure)) {
  print(label_exposure[i])
  list_exposure[[i]] <- separate(list_exposure[[i]], V16, into = c("ID", "REF", "ALT", "rsid", "POS19POS38"), sep = "\t")
  print("column separated")
  colnames(list_exposure[[i]]) <- c("CHR", "POS", "SNPID", "other_allele.exposure", "effect_allele.exposure",
                                    "eaf.exposure", "INFO", "samplesize.exposure", "TEST", 
                                    "beta.exposure", "se.exposure", "CHISQ", "pval.exposure",
                                    "EXTRA", "exposure", "ID", "REF", "ALT", "SNP", "POS19POS38")
  print("colnames changed")
  list_exposure[[i]]$id.exposure <- label_exposure[i]
  print("id.exposure col changed")
  }
## save
exposure <- bind_rows(list_exposure)
write.table(exposure, paste0(DIRECTORY_OUT, "data_exposure.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print("data_exposure.txt saved")

#exposure <- fread(paste0(DIRECTORY_OUT, "data_exposure.txt"))
#list_exposure <- split(exposure, exposure$id.exposure)
#label_exposure <- exposure$id.exposure

# data_outcome ====
list_files_outcome <- list.files("/data/GWAS_data/work/omara_2018_PMID30093612", pattern = "ECAC", full.names = T)
list_files_outcome <- list_files_outcome[3] # noUKBB
## label_outcome
label_outcome <- sub("/data/GWAS_data/work/omara_2018_PMID30093612/", "", list_files_outcome)
label_outcome <- sub(".txt.gz", "", label_outcome)
## extract exposure SNPs from outcome
list_outcome <- list()
print("starting: for loop format outcome")
for (i in 1:length(list_exposure)){
  list_outcome[i] <- lapply(list_files_outcome, 
                            read_outcome_data,
                            snps = list_exposure[[i]]$SNP,
                            sep = "\t",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            se_col = "SE",
                            eaf_col = "EAF",
                            effect_allele_col = "EA",
                            other_allele_col = "OA",
                            pval_col = "P",
                            min_pval = 1e-200,
                            log_pval = FALSE,
                            chr_col = "CHR",
                            pos_col = "POS",
                            phenotype_col = "phenotype")
  print("column changed")
  list_outcome[[i]]$outcome <- label_outcome
  list_outcome[[i]]$id.outcome <- paste0(list_exposure[[i]]$id.exposure[1], "_", list_outcome[[i]]$outcome)
  print("id.outcome col changed")
  }
## save
outcome <- bind_rows(list_outcome)
write.table(outcome, paste0(DIRECTORY_OUT, "data_outcome.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print("data_outcome.txt saved")


# harmonise ====
print("harmonise data")
harmonise_data <- harmonise_data(exposure, outcome, action = 2)
## remove duplicates
harmonise_data$remove_duplicates <- paste0(harmonise_data$SNP, "_", harmonise_data$id.exposure)
harmonise_data <- harmonise_data[!duplicated(harmonise_data$remove_duplicates),]
print("duplicates removed")
## make a list of each harmonised data frame
list_harmonise <- split(harmonise_data, harmonise_data$id.exposure)
## save
write.table(harmonise_data, paste0(DIRECTORY_OUT, "data_harmonised.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print("data_harmonised.txt saved")


## colocalization ====
if ("package:functions" %in% search()) {
  detach("package:functions", unload = TRUE)
}
table_master <- data.frame() # make empty dataframe for final results
lderror <- list()

# start loop for colocalization ====
print("starting: for loop colocalization")
for (i in 1:length(list_harmonise)){
  # make label ====
  label <- paste0(list_harmonise[[i]]$id.exposure[1], ";", list_harmonise[[i]]$outcome[1])
  
  ## make LD matrix ====
  # sometimes no variant extracted from ref file, if error save id and go to next protein
  lderror[[i]] <- try(ld <- ld_matrix_local(
    list_harmonise[[i]]$SNP,
    with_alleles = FALSE, 
    bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
    plink_bin = get_plink_exe())) 
  
  if(class(lderror[[i]]) == "try-error") { 
    print("lderror") 
    next
  }
  
  ld <- ld_matrix_local(
    list_harmonise[[i]]$SNP,
    with_alleles = FALSE, 
    bfile = "/data/GWAS_data/work/references/1kG_v3/EUR/EUR",
    plink_bin = get_plink_exe())
  print("ld matrix made")
  
  # format LD matrix and harmonised list ====
  ld <- ld[which(rownames(ld) %in% list_harmonise[[i]]$SNP), which(colnames(ld) %in% list_harmonise[[i]]$SNP)]
  list_harmonise[[i]] <- list_harmonise[[i]][which(list_harmonise[[i]]$SNP %in% rownames(ld)),]
  ld <- ld[match(list_harmonise[[i]]$SNP,rownames(ld)),]
  ld <- ld[,match(list_harmonise[[i]]$SNP, colnames(ld))]
  list_harmonise[[i]] <- list_harmonise[[i]][match(rownames(ld), list_harmonise[[i]]$SNP),]
  print("ld matrix and harmonised list formatted")
  
  # make lists for coloc ====
  coloc_data_exposure <- list(beta = list_harmonise[[i]]$beta.exposure, 
                              varbeta = list_harmonise[[i]]$se.exposure^2, 
                              MAF = list_harmonise[[i]]$eaf.exposure, 
                              type = "quant", 
                              N = 34557, 
                              snp = rownames(ld), 
                              LD = ld, 
                              position = list_harmonise[[i]]$POS)
  coloc_data_outcome <- list(beta = list_harmonise[[i]]$beta.outcome, 
                             varbeta = list_harmonise[[i]]$se.outcome^2, 
                             MAF = list_harmonise[[i]]$eaf.outcome, 
                             type = "cc", 
                             N = 54884, 
                             s = 0.15957,
                             snp = rownames(ld), 
                             LD = ld, 
                             position = list_harmonise[[i]]$POS)
  
  ## coloc ====  
  # https://chr1swallace.shinyapps.io/coloc-priors/
  # N SNPs = 3000
  print("running coloc.abf")
  coloc_results <- coloc.abf(dataset1 = coloc_data_exposure, 
                             dataset2 = coloc_data_outcome, 
                             p1 = 1E-4, p2 = 1E-4, p12 = 5E-6)
  print("saving pdf")
  pdf(paste0(DIRECTORY_OUT, "figures/", label, ".pdf"), 
      height = 10, width = 10)
  coloc_sensitivity <- my_sensitivity(coloc_results, "H4 > 0.8", 
                                      trait1_title = list_harmonise[[i]]$id.exposure[1], trait2_title = list_harmonise[[i]]$outcome[1])
  dev.off()
  
  # make table ====
  table <- data.frame(
    exposure = list_harmonise[[i]]$exposure[1],
    outcome = list_harmonise[[i]]$outcome[1],
    id = label,
    nsnps = coloc_results["summary"][[1]][1],
    h0 = coloc_results["summary"][[1]][2],
    h1 = coloc_results["summary"][[1]][3],
    h2 = coloc_results["summary"][[1]][4],
    h3 = coloc_results["summary"][[1]][5],
    h4 = coloc_results["summary"][[1]][6])
  
  write.table(table, paste0(DIRECTORY_OUT, "tables/", label, ".txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  table_master <- rbind(table_master, table)
  print("table_master rbinded")
  }
# save ====
list_files_table <- list.files(paste0(DIRECTORY_OUT, "tables/"), full.names = T)
list_table <- lapply(list_files_table, fread, sep = " ", header = F)
for (i in seq_along(list_table)) {
  list_table[[i]] <- separate(list_table[[i]], V1, into = c("exposure", "outcome", "id", "nsnps", "h0", "h1", "h2", "h3", "h4"), sep = "\t")
  list_table[[i]] <- list_table[[i]][-1, ]
  list_table[[i]][ ,4:9] <- as.numeric(list_table[[i]][ ,4:9])
}
table_master <- bind_rows(list_table)
write.table(table_master, paste0(DIRECTORY_OUT, "table_master.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print("table_master saved")