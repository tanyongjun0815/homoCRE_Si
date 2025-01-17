###########################################################################
# Generate overlap matrix (CRE vs Chip-seq peaks)
#          Peaks of ChIP-seq were download from ENCODE.
#          A 0/1 matrix denote overlap between CRE and peaks.
#                                             Tan,Yongjun 20211214,20221227
###########################################################################


##################################### 2022-12-27 #############################################
setwd("/scratch/2022-12-25/bioTanyj/UMI_silencer/R6_sequence_features/CRE_vs_Peaks")
source("/work/bio-tanyj/soft/TanYongjun_code.R")

PeakPath <- "/scratch/2022-12-25/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize" # peaks resize to 1bp
CRE_file <- "../../data_file/SiMerge_EnSTARR_EpiMap_SiReported_Pro_Random_Trans_Kmeans_Cons_20221219.rds"
MinOverlap <- 1 # all ChIP-seq peaks resize to 1bp, means center of peak hit CRE.


## calculate
    # CRE
    cres <- readRDS(CRE_file) %>%
        dplyr::filter(str_detect(ID, "Huang|Pang|Jay|Enhancer_HCT116_Hou", negate = T)) %>%
    dplyr::select(seqnames:end, ID, CRE, cell_line)
    count(cres, CRE, cell_line, ID)

    cres_ls <- cres %>%
        split(., .$cell_line)
    cresGR <- map(cres_ls, GenomicRanges::makeGRangesFromDataFrame) %>%
      map(., GenomicRanges::resize, width = 600, fix = "center")

    library(furrr)
    plan(multisession, workers = 38)

    for(i in names(cres_ls)){
        cat(">>> ", i, "\n")
        peak_ls <- list.files(path = PeakPath, pattern = str_c(i, ".+"))
        peak_ls <- peak_ls[str_detect(peak_ls, "_H2|_H3|_H4|_DHS", negate = T)] # discard epigenetic regions.

        m_tmp <- cres_ls[[i]]

        # # for-loop edition
        # for(j in peak_ls){
        #     cat(j, "; ")
        #     peakGR <- fread(str_c(PeakPath, "/", j), header = F) %>% 
        #         GenomicRanges::makeGRangesFromDataFrame(., seqnames.field = "V1", start.field = "V2", end.field = "V3")
        #     ol <- GenomicRanges::findOverlaps(cresGR[[i]], peakGR, minoverlap = MinOverlap)
        #     m_tmp$tmp <- 0
        #     m_tmp$tmp[unique(ol@from)] <- 1
        #     colnames(m_tmp) <- str_replace_all(colnames(m_tmp), "tmp", str_replace_all(j, ".+?_(.+).bed.gz", "\\1"))
        # }
        # cres_ls[[i]] <- cbind(cres_ls[[i]], m_tmp)


        # purrr edition
        colum_tmp <- rep(0, nrow(cres_ls[[i]]))
        m_tmp <- future_map_dfc(peak_ls, function(x){
            peakGR <- fread(str_c(PeakPath, "/", x), header = F) %>% 
                GenomicRanges::makeGRangesFromDataFrame(., seqnames.field = "V1", start.field = "V2", end.field = "V3")
            ol <- GenomicRanges::findOverlaps(cresGR[[i]], peakGR, minoverlap = MinOverlap)
            tmp <- colum_tmp
            tmp[unique(ol@from)] <- 1
            return(tmp)
        }, .progress = T)

        m_tmp <- as.data.frame(m_tmp)
        colnames(m_tmp) <- str_replace_all(peak_ls, "^.+_(.+)_center.+$", "\\1")
        cres_ls[[i]] <- cbind(cres_ls[[i]], m_tmp)
    }
    plan(sequential)

    saveRDS(cres_ls, file = "CRE_vs_Peaks_Hit_matrix_20221227.rds")


##################################### 1214 #############################################
# setwd("/scratch/2022-02-25/bioTanyj/UMI_silencer/R6_sequence_features/CRE_vs_Peaks")
# source("/work/bio-tanyj/soft/TanYongjun_code.R")
# # PeakPath <- "/scratch/2021-12-07/bioTanyj/silencer_old/S4_TFs/gat/bed_files" # raw chip-seq peaks download from ENCODE/GEO.
# PeakPath <- "/scratch/2022-02-25/bioTanyj/UMI_silencer/R2_ChrState/overlap_GAT/TF_peaks_resize"
# CRE_file <- "../../data_file/SiMerge_EnhancerENCODE_SiReported_Pro_Random_KmeansCluster7_TransCluster4_ConsCluster_20211025.rds"
# # MinOverlap <- 25 # for raw ChIP-seq peaks. The width of some peaks were just about 50bp. (homo_cis/database/ENCODE/Quantiles_of_region_width_Chip-seq_datasets_20201103.csv.gz).
# MinOverlap <- 1 # for all ChIP-seq peaks resize to 10bp.


# ## calculate
#     # CRE
#     cres <- readRDS(CRE_file) %>%
#         dplyr::filter(str_detect(ID, "Huang|Pang|Jay|Enhancer_HCT116_Hou", negate = T)) %>%
#     dplyr::select(seqnames:end, ID, CRE, cell_line)

#     cres_ls <- cres %>%
#         split(., .$cell_line)
#     cresGR <- map(cres_ls, GenomicRanges::makeGRangesFromDataFrame)

#     for(i in names(cres_ls)){
#         cat(">>> ", i, "\n")
#         peak_ls <- list.files(path = PeakPath, pattern = str_c(i, ".+"))

#         m_tmp <- cres_ls[[i]]
#         for(j in peak_ls){
#             cat(j, "; ")
#             peakGR <- fread(str_c(PeakPath, "/", j), header = F) %>% 
#                 GenomicRanges::makeGRangesFromDataFrame(., seqnames.field = "V1", start.field = "V2", end.field = "V3")
#             ol <- GenomicRanges::findOverlaps(cresGR[[i]], peakGR, minoverlap = MinOverlap)
#             m_tmp$tmp <- 0
#             m_tmp$tmp[unique(ol@from)] <- 1
#             colnames(m_tmp) <- str_replace_all(colnames(m_tmp), "tmp", str_replace_all(j, ".+?_(.+).bed.gz", "\\1"))
#         }
#         cres_ls[[i]] <- m_tmp
#     }

#     saveRDS(cres_ls, file = "CRE_vs_Peaks_Hit_matrix_20220228.rds")