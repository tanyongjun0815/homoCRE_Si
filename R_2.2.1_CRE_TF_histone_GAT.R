setwd("/scratch/2022-09-21/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/histone_TFs_resize")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
Num_sample <- "10000"
Num_threads <- 40


###------------Common and cell-specific CRE ------------------------------------------------
# classify enhancers and silencers into two types: Common or cell-specific

# ##  run with gat

# ## collect results. -------------------------
# df <- NULL
# for(i in list.files(pattern = "^GAT.+tsv$")){
#   df <- data.table::fread(i) %>%
#     rbind(., df)
# }
# dim(df)
# data.table::fwrite(df,
#   file = "../GAT_TFs_overlaps_results_all_resizeTFpeaks_10bp_20211213-1.csv.gz")

# ## run GAT of histones -------------------
# peak_path <- "/scratch/2021-12-07/bioTanyj/silencer_old/S4_TFs/gat/bed_files/"

# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".+.bed")
# cre_ls <- cre_ls[str_detect(cre_ls, "CellSpecific|Common")]

# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line <- str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))
#   peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)")]

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", "_20211105"), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, ".tsv"))){
#       cat(">>>>> ", id, "  :::: ", date(), "\n")
#       write.table(NA, file = str_c(id, ".running"))
#       system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
#                   " --segment-file=../bed_files/", i, " ",
#                   " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
#                   " --annotation-file=", peak_path, j, " ",
#                   " -t ", Num_threads, " --random-seed=123456 ",
#                   " --num-samples=", Num_sample,
#                   " --with-segment-tracks ",
#                   " --log=all_logs.log ",
#                   " > ", id, ".tsv", sep = ""))
#     }else{
#     #   cat("***** ", i, "exist", "\n")
#     }

#     # rm lock file
#     if(file.exists(str_c(id, ".running")) & 
#        file.exists(str_c(id, ".tsv"))){
#          system(str_c("rm ", id, ".running", sep = ""))
#        }
#   } # j

# } # i


## run GAT of resized TFs -------------------
peak_path <- "/scratch/2022-09-21/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
cre_ls <- list.files(path = "../bed_files/",
                    pattern = ".+EpiMap.+.bed")

for(i in cre_ls){
  # determine segway file based on cell line
  cell_line <- str_extract(i, "HCT116|K562")

  # list of peak files.
  switch(cell_line,
        "HCT116" = peak_ls <- list.files(path = peak_path,
                                         pattern = str_c(cell_line, ".+", sep = "")),
        "K562" = peak_ls <-  list.files(path = peak_path,
                                         pattern = str_c(cell_line, ".+", sep = "")))

  # discard epigenetic-peaks for KmeansCluster
  if (str_detect(i, "Kmeans")) {
     peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
  }
  id <- str_replace_all(i, ".bed", "")

  # run GAT
  for(j in peak_ls){
    id = str_c("GAT_",
              str_replace_all(j, ".bed.gz", "_"), 
              str_replace_all(i, ".bed", "_20211105"), sep = "")
    
    # run GAT
    if(!file.exists(str_c(id, ".running")) & 
       !file.exists(str_c(id, ".tsv"))){
      cat(">>>>> ", id, "  :::: ", date(), "\n")
      write.table(NA, file = str_c(id, ".running"))
      system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
                  " --segment-file=../bed_files/", i, " ",
                  " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
                  " --annotation-file=", peak_path, j, " ",
                  " -t ", Num_threads, " --random-seed=123456 ",
                  " --num-samples=", Num_sample,
                  " --with-segment-tracks ",
                  " --log=all_logs.log ",
                  " > ", id, ".tsv", sep = ""))
    }else{
    #   cat("***** ", i, "exist", "\n")
    }

    # rm lock file
    if(file.exists(str_c(id, ".running")) & 
       file.exists(str_c(id, ".tsv"))){
         system(str_c("rm ", id, ".running", sep = ""))
       }
  } # j

} # i


# ## chech number of result files. -----------------------------
# peak_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")
# Num_not_run <- 0
# Num_not_finished <- 0
# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", ""), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))){
#       cat(">>>>> ", id, " not Run!\n")
#       Num_not_run <- Num_not_run + 1

#     }else if (file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))) {
#        cat("***** ", id, " not finished!\n")
#       Num_not_finished <- Num_not_finished = 1
#        system(command = str_c("rm ", id, ".running"))
#     }

#   } # j

# } # i

# Num_not_run

# Num_not_finished



##########################################################################################################
###-------------------------------------------DHS region ------------------------------------------------
##### 2021.11.18. Overlap between CREs and TF peaks located in DHS regions.

# # ## run GAT-------------------
# peak_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")

# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", "_20211105"), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, ".tsv"))){
#       cat(">>>>> ", id, "  :::: ", date(), "\n")
#       write.table(NA, file = str_c(id, ".running"))
#       system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
#                   " --segment-file=../bed_files/", i, " ",
#                   " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
#                   " --annotation-file=", peak_path, j, " ",
#                   " -t ", Num_threads, " --random-seed=123456 ",
#                   " --num-samples=", Num_sample,
#                   " --with-segment-tracks ",
#                   " --log=all_logs.log ",
#                   " > ", id, ".tsv", sep = ""))
#     }else{
#       cat("***** ", i, "exist", "\n")
#     }

#     # rm lock file
#     if(file.exists(str_c(id, ".running")) & 
#        file.exists(str_c(id, ".tsv"))){
#          system(str_c("rm ", id, ".running", sep = ""))
#        }
#   } # j

# } # i


# ## chech number of result files. -----------------------------
# peak_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")
# Num_not_run <- 0
# Num_not_finished <- 0
# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", ""), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))){
#       cat(">>>>> ", id, " not Run!\n")
#       Num_not_run <- Num_not_run + 1

#     }else if (file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))) {
#        cat("***** ", id, " not finished!\n")
#       Num_not_finished <- Num_not_finished = 1
#        system(command = str_c("rm ", id, ".running"))
#     }

#   } # j

# } # i

# Num_not_run

# Num_not_finished


##########################################################################################################
###------------resize ChIP-seq peaks------------------------------------------------
# 2021.11.05, peaks of TFs were resize to 10bp.

# ## collect results. -------------------------
# df <- NULL
# for(i in list.files(pattern = "^GAT.+tsv$")){
#   df <- data.table::fread(i) %>%
#     rbind(., df)
# }
# dim(df)
# data.table::fwrite(df, 
#   file = "../GAT_histone_TFs_overlaps_results_all_resizeTFpeaks_10bp_20211105.csv.gz")

# ## resize peaks ----------------------------
# target_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# peak_path <- "/scratch/2021-11-03/bioTanyj/silencer_old/S4_TFs/gat/bed_files/"
# peak_ls <- list.files(path = peak_path, pattern = ".+bed.gz")
# peak_TF <- peak_ls[str_detect(peak_ls, "H2A|H3K|H4|DHS", negate = T)]
# summaryWidth <- NULL
# for(i in peak_TF){
#   tmp <- fread(str_c(peak_path, i, sep = ""), header = F) %>%
#     GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3")
#   summaryWidth <- rbind(summaryWidth, c(i, quantile(width(tmp))))
#   tmp %>%
#     GenomicRanges::resize(width = 10, fix = "center") %>%
#     as.data.frame() %>%
#     dplyr::select(1:3) %>%
#     fwrite(., file = str_c(target_path, str_replace_all(i, ".bed.gz", "_resize.bed.gz")),
#            sep = "\t", col.names = F)
# }
# fwrite(summaryWidth, file = str_c(target_path, "region_width_summary.csv"))


# # ## run GAT-------------------
# peak_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")

# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", "_20211105"), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, ".tsv"))){
#       cat(">>>>> ", id, "  :::: ", date(), "\n")
#       write.table(NA, file = str_c(id, ".running"))
#       system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
#                   " --segment-file=../bed_files/", i, " ",
#                   " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
#                   " --annotation-file=", peak_path, j, " ",
#                   " -t ", Num_threads, " --random-seed=123456 ",
#                   " --num-samples=", Num_sample,
#                   " --with-segment-tracks ",
#                   " --log=all_logs.log ",
#                   " > ", id, ".tsv", sep = ""))
#     }else{
#       cat("***** ", i, "exist", "\n")
#     }

#     # rm lock file
#     if(file.exists(str_c(id, ".running")) & 
#        file.exists(str_c(id, ".tsv"))){
#          system(str_c("rm ", id, ".running", sep = ""))
#        }
#   } # j

# } # i


# ## chech number of result files. -----------------------------
# peak_path <- "/scratch/2021-11-03/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")
# Num_not_run <- 0
# Num_not_finished <- 0
# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", ""), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))){
#       cat(">>>>> ", id, " not Run!\n")
#       Num_not_run <- Num_not_run + 1

#     }else if (file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211105.tsv"))) {
#        cat("***** ", id, " not finished!\n")
#       Num_not_finished <- Num_not_finished = 1
#        system(command = str_c("rm ", id, ".running"))
#     }

#   } # j

# } # i

# Num_not_run

# Num_not_finished


##########################################################################################################
###----------------origin ChIP-seq peaks------------------------------------------------
# 202110, peaks were directly download from ENCODE. NOT used.
# ## collect results. -------------------------
# df <- NULL
# for(i in list.files(pattern = "^GAT.+tsv$")){
#   df <- data.table::fread(i) %>%
#     rbind(., df)
# }
# dim(df)
# data.table::fwrite(df, 
#   file = "../GAT_histone_TFs_overlaps_results_all_20211028.csv.gz")


# # CREs over TFs and histone
# peak_path <- "/scratch/2021-10-25/bioTanyj/silencer_old/S4_TFs/gat/bed_files/"
# cre_ls <- list.files(path = "../bed_files/",
#                     pattern = ".bed")

# ## run GAT-------------------
# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", ""), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211022.tsv"))){
#       cat(">>>>> ", id, "  :::: ", date(), "\n")
#       write.table(NA, file = str_c(id, ".running"))
#       system(paste("/work/bio-tanyj/miniconda3/envs/peaks_overlap/bin/gat-run.py ",
#                   " --segment-file=../bed_files/", i, " ",
#                   " --workspace-file=/data/bio-tanyj/hg19/hg19_chr_size_no_chrXYM_noBlackList_for_GAT.tsv ",
#                   " --annotation-file=", peak_path, j, " ",
#                   " -t ", Num_threads, " --random-seed=123456 ",
#                   " --num-samples=", Num_sample,
#                   " --with-segment-tracks ",
#                   " --log=all_logs.log ",
#                   " > ", id, "_20211022.tsv", sep = ""))
#     }else{
#       cat("***** ", i, "exist", "\n")
#     }

#     # rm lock file
#     if(file.exists(str_c(id, ".running")) & 
#        file.exists(str_c(id, "_20211022.tsv"))){
#          system(str_c("rm ", id, ".running", sep = ""))
#        }
#   } # j

# } # i


# ## chech number of result files. -------------------------

# Num_not_run <- 0
# Num_not_finished <- 0
# for(i in cre_ls){
#   # determine segway file based on cell line
#   cell_line = str_extract(i, "HCT116|K562")

#   # list of peak files.
#   switch(cell_line,
#         "HCT116" = peak_ls <- list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")),
#         "K562" = peak_ls <-  list.files(path = peak_path,
#                                          pattern = str_c(cell_line, ".+", sep = "")))

#   # discard epigenetic-peaks for KmeansCluster
#   if (str_detect(i, "Kmeans")) {
#      peak_ls <- peak_ls[str_detect(peak_ls, "(H2A|H3|H4|DHS)", negate = T)]
#   }
#   id <- str_replace_all(i, ".bed", "")

#   # run GAT
#   for(j in peak_ls){
#     id = str_c("GAT_",
#               str_replace_all(j, ".bed.gz", "_"), 
#               str_replace_all(i, ".bed", ""), sep = "")
    
#     # run GAT
#     if(!file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211022.tsv"))){
#       cat(">>>>> ", id, " not Run!\n")
#       Num_not_run <- Num_not_run + 1

#     }else if (file.exists(str_c(id, ".running")) & 
#        !file.exists(str_c(id, "_20211022.tsv"))) {
#        cat("***** ", id, " not finished!\n")
#       Num_not_finished <- Num_not_finished = 1
#       #  system(command = str_c("rm ", id, ".running"))
#     }

#   } # j

# } # i

# Num_not_run

# Num_not_finished


