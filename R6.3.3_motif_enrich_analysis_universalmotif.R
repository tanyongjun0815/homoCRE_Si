#############################################################################
# Identify motifs enriched in CREs use other CREs as background.
#                       Tan, Yongjun 20211203, 20221227, 20230412
#############################################################################
  rm(list = ls())
  source("/work/bio-tanyj/soft/TanYongjun_code.R")
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(universalmotif)
  library(tidyverse)
  setwd("/scratch/2023-04-18/bioTanyj/UMI_silencer/R6_sequence_features/motif_enrichment_analysis/universalmotif_package")
  NumThreads <- 15 # for furrr.

##################################################### 2022-12-27 #########################################################
# prepare DNAStringSet (in a list) of CREs
  cres <- read_rds("/scratch/2023-04-18/bioTanyj/UMI_silencer/data_file/SiMerge_EnSTARR_EpiMap_SiReported_Pro_Random_Trans_Kmeans_Cons_20221219.rds") %>%
    dplyr::filter(str_detect(ID, "Huang|Jaya|Pang|Enhancer_HCT116_Hou", negate = T)) %>%
    left_join(., fread("/scratch/2023-04-18/bioTanyj/UMI_silencer/data_file/Append_info_Strand_cellline_TSS_ProHit_20221218.csv.gz")) # nolint
  count(cres, ID)

  # each CRE as a whole.
  timestamp()
  cat(">>> truncate DNA sequence.\n")
  cresFA_ls <- cres %>%
    dplyr::select(-strand) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  # Kmeans cluster
  KmeansFA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random")) %>%
    dplyr::mutate(ID = str_c(ID, "_Kmeans", KmeansCluster)) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  # Kmeans cluster C136 vs C2457
  KmeansGroupFA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random")) %>%
    dplyr::mutate(ID = ifelse(str_detect(KmeansCluster, "C1|C3|C6"),
                              str_c(ID, "_KmeansC136"),
                              str_c(ID, "_KmeansC457"))) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  # Trans group
  TransFA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random")) %>%
    dplyr::mutate(ID = str_c(ID, "_Trans", as.numeric(TransGroup_EpiMap))) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)
  
  # Trans group merg (T3+T4+T5 and T4+T5)
  Trans345FA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random")) %>%
    dplyr::mutate(ID = str_c(ID, "_Trans", as.numeric(TransGroup_EpiMap))) %>%
    dplyr::mutate(ID = str_replace_all(ID, c("Trans3$" = "Trans345",
                                              "Trans4$" = "Trans345",
                                              "Trans5$" = "Trans345"))) %>%
    dplyr::filter(str_detect(ID, "Trans345")) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  Trans45FA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random")) %>%
    dplyr::mutate(ID = str_c(ID, "_Trans", as.numeric(TransGroup_EpiMap))) %>%
    dplyr::mutate(ID = str_replace_all(ID, c("Trans4$" = "Trans45",
                                              "Trans5$" = "Trans45"))) %>%
    dplyr::filter(str_detect(ID, "Trans45")) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  # Kmeans cluster X Trans group (Only high and low trans group)
  TransKmeansFA_ls <- cres %>%
    dplyr::select(-strand) %>%
    dplyr::mutate(TransGroup_EpiMap = str_c("TransG", as.numeric(TransGroup_EpiMap))) %>%
    dplyr::filter(str_detect(ID, "Silencer|Enhancer|Random"),
                  str_detect(ID, "EpiMap", negate = T),
                  str_detect(TransGroup_EpiMap, "TransG2|TransG3", negate = T)) %>% # count(ID, cell_line, TransGroup_EpiMap, KmeansCluster) %>% pivot_wider(names_from = KmeansCluster, values_from = n) %>% kable
    dplyr::mutate(TransGroup_EpiMap = str_replace_all(TransGroup_EpiMap, "TransG4|TransG5", "TransG45"), # nolint
                  ID = str_c("Intersection", ID, "K", KmeansCluster,
                  TransGroup_EpiMap, sep = "_")) %>%
    split(., .$ID) %>%
    map(., makeGRangesFromDataFrame) %>%
    map(., GenomicRanges::resize, width = 600, fix = "center") %>%
    map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

  # discard sequences contain too many N.
  DisCardSeqWithNs <- function(FA_ls){
    NuclFreqLs <- map(FA_ls, alphabetFrequency, as.prob = T) %>%
      map(., function(x){
        which(x[,which(colnames(x) == "N")] < 0.95)
      })
    for (i in names(FA_ls)) {
      FA_ls[[i]] <- FA_ls[[i]][NuclFreqLs[[i]]]
    }
    return(FA_ls)
  }
  cresFA_ls <- DisCardSeqWithNs(cresFA_ls)
  clusterFA_ls <- DisCardSeqWithNs(c(KmeansFA_ls, KmeansGroupFA_ls, TransFA_ls,
                                     Trans345FA_ls, Trans45FA_ls, TransKmeansFA_ls))

  # rm(KmeansFA_ls, TransFA_ls)
  map_int(clusterFA_ls, length)

# motifs
    # PWM_matrix_all <- read_rds("/scratch/2022-03-01/bioTanyj/UMI_silencer/R6_sequence_features/Affinity_score_matrix/PWMs_matrix_5594_motifs_all_20211130.rds")
    PWM_S4_all <- read_rds("/scratch/2023-04-18/bioTanyj/UMI_silencer/R6_sequence_features/Affinity_score_matrix/PWMs_S4_PWMEnrich_1284_best_motifs_20221227.rds")

# combination of CRE pairs
    cresAll_FA_all <- c(cresFA_ls, clusterFA_ls)
    # rm(cresFA_all, clusterFA_ls)
    names(cresAll_FA_all)

    # (Si vs En; HCT116 vs K562)
    df <- rbind(c("Enhancer_HCT116_STARR-seq", "Silencer_HCT116_Hou"),
                c("Enhancer_K562_STARR-seq", "Silencer_K562_Hou"),
                c("Enhancer_HCT116_EpiMap", "Silencer_HCT116_Hou"),
                c("Enhancer_K562_EpiMap", "Silencer_K562_Hou"),
                c("Enhancer_HCT116_STARR-seq", "Enhancer_K562_STARR-seq"),
                c("Silencer_HCT116_Hou", "Silencer_K562_Hou")) %>%
                as.data.frame() %>%
          dplyr::mutate(V1 = as.character(V1), V2 = as.character(V2))
    df <- data.frame(V1 = c(as.character(df$V1), as.character(df$V2)), V2 = c(as.character(df$V2), as.character(df$V1)))

    # CRE vs Random
    df <- data.frame(V1 = c(unique(c(df$V1, df$V2)), "Promoter_HCT116_UCSC", "Promoter_K562_UCSC"),
                           V2 = NA) %>%
                dplyr::mutate(V2 = ifelse(str_detect(V1, "HCT116"), "Random_HCT116_Tan", "Random_K562_Tan")) %>%
                rbind(., df)

    # KmeansCluster/KmeansClusterGroup/Trans/Intersection Grouop vs Random (the corresponding group)
    group_ls <- names(clusterFA_ls)
    group_ls[str_detect(group_ls, "Random")]
    group_ls <- group_ls[str_detect(group_ls, "Kmeans|Trans")]
    group_ls <- group_ls[str_detect(group_ls, "Random", negate = T)]

    df <- data.frame(V1 = group_ls, 
                     tmp = str_replace_all(group_ls, "^.+(Hou|seq|EpiMap)_(.+?)$", "\\2")) %>%
          dplyr::mutate(V2 = ifelse(str_detect(V1, "HCT116"), str_c("Random_HCT116_Tan_", tmp), str_c("Random_K562_Tan_", tmp)),
                        V2 = ifelse(str_detect(V1, "Intersection_"), str_c("Intersection_", V2), V2)) %>%
          dplyr::select(-tmp) %>%
          rbind(., df)

    # KmeansCluster/Trans/Intersection Grouop vs Random (All, 20230412)
    df <- data.frame(V1 = group_ls) %>%
        dplyr::mutate(V2 = ifelse(str_detect(V1, "HCT116"), "Random_HCT116_Tan", "Random_K562_Tan")) %>%
        rbind(., df)

    # Random vs Random (estimate false positve rate)
    df <- rbind(df,
                data.frame(V1 = c("Random_HCT116_Tan", "Random_K562_Tan"),
                           V2 = c("Random_K562_Tan", "Random_HCT116_Tan")))

    fwrite(df, file = "CRE_combination_list.csv")

    # run
    timestamp()
    cat(">>> Run enrichment analysis.\n")
    library(furrr)
    plan(multisession, workers = NumThreads)
    for (i in 1:nrow(df)) {
      id1 <- df$V1[i]
      id2 <- df$V2[i]
      lockID <- str_c(id1, id2, ".lock", sep = "_")
      resultID <- str_c(id1, id2, "universalmotif_enrich_motif.rds", sep = "_")

      if (!file.exists(lockID) & !file.exists(resultID)) {
        write(NA, file = lockID)
        timestamp()
        cat(">>> ", lockID, "\n")
        seqTMP <- cresAll_FA_all[[id1]]
        bkg.seqTMP <- cresAll_FA_all[[id2]]
        tmp <- future_map_dfr(PWM_S4_all,
            function(x){
              enrich_motifs(motifs = x,
                            sequences = seqTMP,
                            bkg.sequences = bkg.seqTMP,
                            RC = T) %>%
                as.data.frame() %>%
                return()},
            .id = "ID")
        saveRDS(tmp, file = resultID)
        system(str_c("rm ", lockID))
      } # if
      
    } # i-for
    plan(sequential)



# ##################################################### 2022-12-27 #########################################################
# # prepare DNAStringSet (in a list) of CREs
#   cres <- read_rds("/scratch/2022-03-01/bioTanyj/UMI_silencer/data_file/SiMerge_EnhancerENCODE_SiReported_Pro_Random_KmeansCluster7_TransCluster4_ConsCluster_20211025.rds") %>%
#     dplyr::filter(str_detect(ID, "Huang|Jaya|Pang|Enhancer_HCT116_Hou", negate = T)) %>%
#     left_join(., fread("/scratch/2022-03-01/bioTanyj/UMI_silencer/data_file/Append_info_strand_cell_specific_20211209.csv.gz"))
  
#   # each CRE as a whole.
#   cresFA_ls <- cres %>%
#     dplyr::select(-strand) %>%
#     split(., .$ID) %>%
#     map(., makeGRangesFromDataFrame) %>%
#     map(., GenomicRanges::resize, width = 600, fix = "center") %>%
#     map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

#   # Kmeans cluster and Common/Cellspecific
#   KmeansFA_ls <- cres %>%
#     dplyr::select(-strand) %>%
#     dplyr::filter(str_detect(ID, "Silencer|Enhancer|Promoter")) %>%
#     dplyr::mutate(ID = str_c(ID, "_Kmeans", Kmeans_cluster)) %>%
#     split(., .$ID) %>%
#     map(., makeGRangesFromDataFrame) %>%
#     map(., GenomicRanges::resize, width = 600, fix = "center") %>%
#     map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

#   CommonSpeciFA_ls <- cres %>%
#     dplyr::select(-strand) %>%
#     dplyr::filter(str_detect(ID, "Silencer|Enhancer")) %>%
#     dplyr::mutate(ID = str_c(ID, "_", cell_line_type)) %>%
#     split(., .$ID) %>%
#     map(., makeGRangesFromDataFrame) %>%
#     map(., GenomicRanges::resize, width = 600, fix = "center") %>%
#     map(., BSgenome::getSeq, x = BSgenome.Hsapiens.UCSC.hg19)

#   # discard sequences contain too many N.
#   DisCardSeqWithNs <- function(FA_ls){
#     NuclFreqLs <- map(FA_ls, alphabetFrequency, as.prob = T) %>%
#       map(., function(x){
#         which(x[,which(colnames(x) == "N")] < 0.95)
#       })
#     for (i in names(FA_ls)) {
#       FA_ls[[i]] <- FA_ls[[i]][NuclFreqLs[[i]]]
#     }
#     return(FA_ls)
#   }
#   cresFA_ls <- DisCardSeqWithNs(cresFA_ls)
#   clusterFA_ls <- DisCardSeqWithNs(c(KmeansFA_ls, CommonSpeciFA_ls))
#   rm(KmeansFA_ls, CommonSpeciFA_ls)

# # motifs
#     # PWM_matrix_all <- read_rds("/scratch/2022-03-01/bioTanyj/UMI_silencer/R6_sequence_features/Affinity_score_matrix/PWMs_matrix_5594_motifs_all_20211130.rds")
#     PWM_S4_all <- read_rds("/scratch/2022-03-01/bioTanyj/UMI_silencer/R6_sequence_features/Affinity_score_matrix/PWMs_S4_PWMEnrich_5594_motifs_all_20211130.rds")


#   # combination 
#     # (Si vs En; HCT116 vs K562)
#     df <- rbind(c("Enhancer_HCT116_ENCODE", "Silencer_HCT116_Hou"),
#                 c("Enhancer_K562_ENCODE", "Silencer_K562_Hou"),
#                 c("Enhancer_HCT116_ENCODE", "Enhancer_K562_ENCODE"),
#                 c("Silencer_HCT116_Hou", "Silencer_K562_Hou")) %>%
#                 as.data.frame() %>%
#           dplyr::mutate(V1 = as.character(V1), V2 = as.character(V2))
#     df <- data.frame(V1 = c(as.character(df$V1), as.character(df$V2)), V2 = c(as.character(df$V2), as.character(df$V1)))

#     # CRE vs Random
#     df <- rbind(df, 
#                 data.frame(V1 = unique(c(df$V1, df$V2)),
#                            V2 = c("Random_HCT116_Tan", "Random_K562_Tan", "Random_HCT116_Tan", "Random_K562_Tan")))

#     # (KmeansCluster/Common/CellSpecific vs Random)
#     df <- rbind(df,
#                 data.frame(V1 = names(clusterFA_ls),
#                            V2 = ifelse(str_detect(names(clusterFA_ls), "HCT116"),
#                                         "Random_HCT116_Tan", "Random_K562_Tan")))

#     # Random vs Random (estimate false positve rate)
#     df <- rbind(df,
#                 data.frame(V1 = c("Random_HCT116_Tan", "Random_K562_Tan"),
#                            V2 = c("Random_K562_Tan", "Random_HCT116_Tan")))

#     cresFA_ls <- c(cresFA_ls, clusterFA_ls)

#     # run
#     library(furrr)
#     plan(multisession, workers = 40)
#     for (i in 1:nrow(df)) {
#       id1 <- df$V1[i]
#       id2 <- df$V2[i]
#       lockID <- str_c(id1, id2, ".lock", sep = "_")
#       resultID <- str_c(id1, id2, "universalmotif_enrich_motif.rds", sep = "_")

#       if (!file.exists(lockID) & !file.exists(resultID)) {
#         write(NA, file = lockID)
#         timestamp()
#         cat(">>> ", lockID, "\n")
#         seqTMP <- cresFA_ls[[id1]]
#         bkg.seqTMP <- cresFA_ls[[id2]]
#         tmp <- future_map_dfr(PWM_S4_all,
#             function(x){
#               enrich_motifs(motifs = x,
#                             sequences = seqTMP,
#                             bkg.sequences = bkg.seqTMP,
#                             RC = T) %>%
#                 as.data.frame() %>%
#                 return()},
#             .id = "ID")
#         saveRDS(tmp, file = resultID)
#         system(str_c("rm ", lockID))
#       } # if
      
#     } # i-for
#     plan(sequential)
