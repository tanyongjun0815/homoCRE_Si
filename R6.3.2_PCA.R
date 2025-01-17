#############################################################################
# Perform PCA based on mean affinity score matrix of CRE in two cell lines
# Analysis:
#     1. PCA using function from package stats and factoextra.
#     2. some seqs were discard in calculate mean affinity score for to much Ns.
#           (need ~ 80GB RAM for each cell line)
#                       Tan, Yongjun 20211204,20230421
#############################################################################

  # source("/work/bio-tanyj/soft/TanYongjun_code.R")
  library(tidyverse)
  library(factoextra)
  setwd("/scratch/2023-04-18/bioTanyj/UMI_silencer/R6_sequence_features/Affinity_score_matrix/means_score_matrix")

  # cres_filt <- readRDS("CRE_list_affinity_score_20211130.rds") 
  cres_filt <- readRDS("CRE_list_with_affinity_score_20221227.rds") 
  unique(cres_filt$ID)
  # score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_20211130.rds") %>%
  score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_20221227.rds") %>%
  # score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_sampled_20221227.rds") %>%
    as.matrix()

  # name of sequences
  row.names(score_matrix) <- cres_filt$name

 # discard reported CREs.
  cres_filt <- cres_filt %>%
    dplyr::filter(str_detect(ID, "HepG2|Huang|Jay|Pang|Random", negate = T))
# ## Only best motifs (328604 seqs X 1496 motifs)-------------------------------------------
#  # PCA of all CREs in two cell lines.
#   MotifLs <- read.csv("../Human_TF_MotifList_v_1.01.csv") %>%
#     dplyr::filter(!is.na(`Best.Motif.s....Figure.2A.`))

#   M_tmp <- score_matrix[, colnames(score_matrix) %in% MotifLs$CIS.BP.ID]

#   # PCA
#   cat(">>> Perform PCA, Number of CREs:", nrow(M_tmp), "\n")
#   timestamp()
#   res <- prcomp(M_tmp, scale. = T, rank. = 6)
  
#   pc_var <- get_pca_var(res)
#   pc_ind <- get_pca_ind(res)
#   res_summary <- summary(res)
  
#   saveRDS(pc_var, file = str_c("PCA_results_HCT116_K562_motif_score_BestMotif_PCA_20220315.rds", sep = ""))
#   saveRDS(pc_ind, file = str_c("PCA_results_HCT116_K562_CREs_score_BestMotif_PCA_20220315.rds", sep = ""))
#   saveRDS(res_summary, file = str_c("PCA_summary_HCT116_K562_PCs_importance_BestMotif_PCA_20220315.rds", sep = ""))


 # PCA in each cell line
  MotifLs <- read.csv("../Human_TF_MotifList_v_1.01.csv") %>%
    dplyr::filter(!is.na(`Best.Motif.s....Figure.2A.`))
  M_best <- score_matrix[, colnames(score_matrix) %in% MotifLs$CIS.BP.ID]


  for(i in c("HCT116", "K562")){
      if (!file.exists(str_c("PCA_", i, "_20230414.calculating", sep = ""))) {
          write.csv(NA, file = str_c("PCA_", i, "_20230414.calculating", sep = ""))
          timestamp()
          cat(">>> Prepare matrix for PCA:")
          cat(">>>>", i, "\n")
          # PCA
          cres_tmp <- cres_filt %>%
            dplyr::filter(cell_line == i)
          M_tmp <- M_best[row.names(M_best) %in% cres_tmp$name,]
          cat(">>> Perform PCA of ", i, "; Number of CREs:", nrow(M_tmp), "\n")
          timestamp()
          res <- prcomp(M_tmp, scale. = T, rank. = 3)
          
          pc_var <- get_pca_var(res)
          pc_ind <- get_pca_ind(res)
          res_summary <- summary(res)
          
          # saveRDS(pc_var, file = str_c("PCA_results_", i, "_motif_score_BestMotif_PCA_20220315.rds", sep = ""))
          # saveRDS(pc_ind, file = str_c("PCA_results_", i, "_CREs_score_BestMotif_PCA_20220315.rds", sep = ""))
          # saveRDS(res_summary, file = str_c("PCA_summary_", i, "_PCs_importance_BestMotif_PCA_20220315.rds", sep = ""))
          saveRDS(pc_var, file = str_c("PCA_results_", i, "_motif_score_BestMotif_PCA_20230412.rds", sep = ""))
          saveRDS(pc_ind, file = str_c("PCA_results_", i, "_CREs_score_BestMotif_PCA_20230412.rds", sep = ""))
          saveRDS(res_summary, file = str_c("PCA_summary_", i, "_PCs_importance_BestMotif_PCA_20230412.rds", sep = ""))
      } # if
  } # for


# ## all (328604 seqs X 5594 motifs)-------------------------------------------

#   M_tmp <- score_matrix
  
#   # PCA
#   cat(">>> Perform PCA, Number of CREs:", nrow(M_tmp), "\n")
#   timestamp()
#   res <- prcomp(M_tmp, scale. = T, rank. = 6)
  
#   pc_var <- get_pca_var(res)
#   pc_ind <- get_pca_ind(res)
#   res_summary <- summary(res)
  
#   saveRDS(pc_var, file = str_c("PCA_results_HCT116_K562_motif_score_AllMotifs_20211205.rds", sep = ""))
#   saveRDS(pc_ind, file = str_c("PCA_results_HCT116_K562_CREs_score_AllMotifs_20211205.rds", sep = ""))
#   saveRDS(res_summary, file = str_c("PCA_summary_HCT116_K562_PCs_importance_AllMotifs_20211205.rds", sep = ""))
  

## Hou-lab C ------------------------------------------------------------------
  source("/data/tanyj/scripts/TanYongjun_code.R")
  library(tidyverse)
  library(factoextra)
  setwd("/data/tanyj/silencer/A6_sequences/motif_hit_matrix/mean_affinity_score_matrix")

  # cres_filt <- readRDS("CRE_list_affinity_score_20211130.rds") 
  cres_filt <- readRDS("CRE_list_with_affinity_score_20221227.rds") 
  unique(cres_filt$ID)
  # score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_20211130.rds") %>%
  score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_20221227.rds") %>%
  # score_matrix <- readRDS("Mean_affinity_score_matrix_motifScore_sampled_20221227.rds") %>%
    as.matrix()

  # name of sequences
  row.names(score_matrix) <- cres_filt$name

 # discard reported CREs.
  cres_filt <- cres_filt %>%
    dplyr::filter(str_detect(ID, "HepG2|Huang|Jay|Pang|Random", negate = T))

 # PCA in each cell line

  for(i in c("HCT116", "K562")){
      if (!file.exists(str_c("PCA_", i, "_20230414.calculating", sep = ""))) {
          write.csv(NA, file = str_c("PCA_", i, "_20230414.calculating", sep = ""))
          timestamp()
          cat(">>> Prepare matrix for PCA:")
          cat(">>>>", i, "\n")
          # PCA
          cres_tmp <- cres_filt %>%
            dplyr::filter(cell_line == i)
          M_tmp <- score_matrix[row.names(score_matrix) %in% cres_tmp$name,]
          cat(">>> Perform PCA of ", i, "; Number of CREs:", nrow(M_tmp), "\n")
          timestamp()
          res <- prcomp(M_tmp, scale. = T, rank. = 3)
          
          pc_var <- get_pca_var(res)
          pc_ind <- get_pca_ind(res)
          res_summary <- summary(res)
          
          # saveRDS(pc_var, file = str_c("PCA_results_", i, "_motif_score_BestMotif_PCA_20220315.rds", sep = ""))
          # saveRDS(pc_ind, file = str_c("PCA_results_", i, "_CREs_score_BestMotif_PCA_20220315.rds", sep = ""))
          # saveRDS(res_summary, file = str_c("PCA_summary_", i, "_PCs_importance_BestMotif_PCA_20220315.rds", sep = ""))
          saveRDS(pc_var, file = str_c("PCA_results_", i, "_motif_score_BestMotif_PCA_20230412.rds", sep = ""))
          saveRDS(pc_ind, file = str_c("PCA_results_", i, "_CREs_score_BestMotif_PCA_20230412.rds", sep = ""))
          saveRDS(res_summary, file = str_c("PCA_summary_", i, "_PCs_importance_BestMotif_PCA_20230412.rds", sep = ""))
      } # if
  } # for
