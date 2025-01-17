##  run with gat
setwd("/scratch/2023-04-28/bioTanyj/UMI_silencer/R7_eQTL_GWAS/GAT")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
Num_sample <- "10000"
Num_threads <- 10

# CREs over Segway(merged states) --------------------
crePath <- "/scratch/2023-04-28/bioTanyj/UMI_silencer/data_file/bed_files"
SNP_file <- "/scratch/2023-04-28/bioTanyj/UMI_silencer/R7_eQTL_GWAS/sig_variants_eQTL_GWAS_hg19_20220205.bed.gz"
cre_ls <- list.files(path = crePath,
                    pattern = ".+600bp.bed")
cre_ls <- cre_ls[str_detect(cre_ls, "Intersection|EpiMap_Kmeans", negate = T)]
cre_ls

Sys.sleep(sample((1:1000/100), 1)) # prevent repeat running GAT for the same bed file.

for(i in cre_ls){

  id <- str_replace_all(i, ".bed", "")
  # run GAT
  if(!file.exists(str_c(id, ".running", sep = "")) &
     !file.exists(str_c("GAT_", id, "20221220.tsv", sep = ""))){
    cat(">>>>> ", id, "  :::: ", date(), "\n")
    write.table(NA, file = str_c(id, ".running"))
    system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
                " --segment-file=", crePath, "/", i, " ",
                " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
                " --annotation-file=", SNP_file, " ",
                " -t ", Num_threads, " --random-seed=123456 ",
                " --num-samples=", Num_sample,
                " --with-segment-tracks ",
                " --log=all_logs.log ",
                " > ", "GAT_", id, "20221220.tsv", sep = ""))

  }else{
    cat("***** ", i, "exist", "\n")
  }

    if(file.exists(str_c(id, ".running", sep = "")) &
     file.exists(str_c("GAT_", id, "20221220.tsv", sep = ""))){
       system(command = str_c("rm ", id, ".running", sep = ""))
     }

}


## collect results. -------------------------
# # merged
df <- NULL
for(i in list.files(pattern = "^GAT.+.tsv$")){
  df <- data.table::fread(i) %>%
    rbind(., df)
}
dim(df)
data.table::fwrite(df, 
  file = "../GAT_CRE_eQTL_GWAS_overlaps_results_all_20221221.csv.gz")
