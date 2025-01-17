##  run with gat
setwd("/scratch/2023-05-29/bioTanyj/UMI_silencer/R4_Conservation/TEs")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
Num_sample <- "10000"
Num_threads <- 5

# silencer vs enhancer table
cre_path <- "../../data_file/bed_files/"
cre_ls <- list.files(path = cre_path,
                    pattern = ".+1bp.bed")
cre_ls <- cre_ls[str_detect(cre_ls, "Intersection", negate = T)]
cre_ls

Sys.sleep(sample((1:1000/100), 1)) # prevent repeat running GAT for the same bed file.

## run GAT-------------------
for(i in cre_ls){

    for(j in list.files(path = "./", pattern = "^Repeat.+bed.gz$")){
        id <- str_c(str_replace_all(j, ".+_(.+)_2021.+", "\\1"),
                    str_replace_all(i, ".bed", ""),
                    sep = "_")
        # run GAT
        if(!file.exists(str_c(id, ".running")) & 
          !file.exists(str_c(id, "_20221220.tsv"))){
          cat(">>>>> ", id, "  :::: ", date(), "\n")
          write.table(NA, file = str_c(id, ".running"))
          system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
                      " --segment-file=", cre_path, i, " ",
                      " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
                      " --annotation-file=./", j, " ",
                      " -t ", Num_threads, " --random-seed=123456 ",
                      " --num-samples=", Num_sample,
                      " --with-segment-tracks ",
                      " --log=all_logs.log ",
                      " > ", id, "_20221220.tsv", sep = ""))
        }else{
          cat("***** ", i, "exist", "\n")
        }

        # rm lock file
        if(file.exists(str_c(id, ".running")) & 
          file.exists(str_c(id, "_20221220.tsv"))){
            system(str_c("rm ", id, ".running", sep = ""))
          }
    } # j

} # i


## collect results. -------------------------
df <- NULL
for(i in list.files(pattern = "^rep.+tsv$")){
  df <- data.table::fread(i) %>%
    dplyr::mutate(ID = str_extract(i, "repName|repFamily")) %>%
    rbind(., df)
}
dim(df)

data.table::fwrite(df, 
  file = "./GAT_CRE_TEs_overlaps_results_all_20230530.csv.gz")