##  run with gat
setwd("/scratch/2022-12-20/bioTanyj/UMI_silencer/R5_Distribution/GAT")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
Num_sample <- "10000"
Num_threads <- 15

# silencer vs enhancer table
cre_path <- "/scratch/2022-12-20/bioTanyj/UMI_silencer/data_file/bed_files/"
cre_ls <- list.files(path = cre_path,
                    pattern = ".+1bp.bed")
cre_ls <- cre_ls[str_detect(cre_ls, "Promoter", negate = T)]
cre_ls
Sys.sleep(sample(1:1000, 1)/ 100)

## run GAT-------------------
for(i in cre_ls){

    id <- str_c(str_replace_all(i, ".bed", ""),
                sep = "_")
    # run GAT
    if(!file.exists(str_c(id, ".running")) & 
      !file.exists(str_c(id, "_20221220.tsv"))){
      cat(">>>>> ", id, "  :::: ", date(), "\n")
      write.table(NA, file = str_c(id, ".running"))
      system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
                  " --segment-file=", cre_path, i, " ",
                  " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
                  " --annotation-file=./Hs_GenomicFeatures_20211105.bed.gz",
                  " -t ", Num_threads, " --random-seed=123456 ",
                  " --num-samples=", Num_sample,
                  " --with-segment-tracks ",
                  " --log=all_logs.log ",
                  " > ", id, "_20221220.tsv", sep = ""))
                  
      # rm lock file
      if(file.exists(str_c(id, ".running")) & 
        file.exists(str_c(id, "_20221220.tsv"))){
          system(str_c("rm ", id, ".running", sep = ""))
        }
    # rm lock file
    if(file.exists(str_c(id, ".running")) & 
      file.exists(str_c(id, "_20221220.tsv"))){
        system(str_c("rm ", id, ".running", sep = ""))
      }
    }else{
      cat("***** ", i, "exist", "\n")
    }

} # i


## collect results. -------------------------
df <- NULL
for(i in list.files(pattern = ".+tsv$")){
  df <- data.table::fread(i) %>%
    dplyr::mutate(ID = str_replace_all(i, "_20221220.tsv", "")) %>%
    rbind(., df)
}
dim(df)

data.table::fwrite(df, 
  file = "./GAT_CRE_GenomicFeatures_overlaps_results_all_20221220.csv.gz")