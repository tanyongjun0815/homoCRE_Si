###################################################################################
# Identify TFs enriched in loop anchors with GAT.
#     four types of loop anchors: hit by Si/Pro/En/Random.
#     all loop anchors in each cell line were used as workspace.
#                                   Tan,Yongjun 20220305
###################################################################################
##  run with gat
setwd("/scratch/2022-03-01/bioTanyj/UMI_silencer/R5_Distribution/TFs_peaks_enriched_loop_anchor")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
Num_sample <- "10000"
Num_threads <- 40


# ## run GAT-------------------
peak_path <- "/scratch/2022-03-01/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
cre_ls <- list.files(path = "./", pattern = ".+(Silencer|Enhancer|Promoter|Silencer|Random).+.bed")

for(i in cre_ls){
  # determine segway file based on cell line
  cell_line = str_extract(i, "HCT116|K562")

  # list of peak files.
  switch(cell_line,
        "HCT116" = workSpace <- "/scratch/2022-03-01/bioTanyj/UMI_silencer/R5_Distribution/TFs_peaks_enriched_loop_anchor/loops_anchor_HCT116_all_10kbWidth.bed",
        "K562" = workSpace <-  "/scratch/2022-03-01/bioTanyj/UMI_silencer/R5_Distribution/TFs_peaks_enriched_loop_anchor/loops_anchor_K562_all_10kbWidth.bed")

  # all peaks anchors as work space in GAT.
  switch(cell_line,
        "HCT116" = peak_ls <- list.files(path = peak_path,
                                         pattern = str_c(cell_line, ".+", sep = "")),
        "K562" = peak_ls <-  list.files(path = peak_path,
                                         pattern = str_c(cell_line, ".+", sep = "")))

  # run GAT
  for(j in peak_ls){
    id = str_c("GAT_",
              str_replace_all(j, "resize.bed.gz", ""), 
              str_replace_all(i, "loops_anchor_(.+)_10kb.+", "\\1"), sep = "")

    # run GAT
    if(!file.exists(str_c(id, ".running")) & 
       !file.exists(str_c(id, ".tsv"))){
      cat(">>>>> ", id, "  :::: ", date(), "\n")
      write.table(NA, file = str_c(id, ".running"))
      system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
                  " --segment-file=./", i, " ",
                  " --workspace-file=", workSpace, 
                  " --annotation-file=", peak_path, j, " ",
                  " -t ", Num_threads, " --random-seed=123456 ",
                  " --num-samples=", Num_sample,
                  " --with-segment-tracks ",
                  " --log=all_logs.log ",
                  " > ", id, ".tsv", sep = ""))
    }else{
      cat("***** ", i, "exist", "\n")
    }

    # rm lock file
    if(file.exists(str_c(id, ".running")) & 
       file.exists(str_c(id, ".tsv"))){
         system(str_c("rm ", id, ".running", sep = ""))
       }
  } # j

} # i



# ## collect results. -------------------------
df <- NULL
for(i in list.files(pattern = "^GAT.+tsv$")){
  df <- data.table::fread(i) %>%
    rbind(., df)
}
dim(df)
data.table::fwrite(df, 
  file = "./GAT_TFs_Anchors_overlaps_results_all_resizeTFpeaks_10bp_2020304.csv.gz")



