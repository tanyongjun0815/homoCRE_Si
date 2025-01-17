    setwd("/scratch/2022-12-14/bioTanyj/UMI_silencer/R3.Transition/EpiMap_Hit_Matrix")
    source("/work/bio-tanyj/scripts/TanYongjun_code.R", echo = F)

    cre_all <- read_rds("./SiMerge_EnhancerENCODE_Epimap_SiReported_Pro_Random_20220803.rds")
    creGR <- makeGRangesFromDataFrame(cre_all) %>%
        GenomicRanges::resize(width = 1, fix = "center")
    # SelectMeta <- fread("/data/bio-tanyj/hg19/EpiMap/Representative_100_biosamples_hiearchical_cluster_20220804.csv.gz") # 100 representative samples with arbitrary number of clusters.
    SelectMeta <- fread("/data/bio-tanyj/hg19/EpiMap/Representative_160_biosamples_hiearchical_cluster_DynamicCut_20221218.csv.gz")
    NumThreds <- 38
    
    # Enhancers matrix
    bed_ls <- list.files(path = "/data/bio-tanyj/hg19/EpiMap/Enhancers/",
                        pattern = "BSS.+bed.gz") %>%
        as.list()
    names(bed_ls) <- str_replace_all(bed_ls, "(^.+?)_.+", "\\1")
    
    HitEn <- function(x){
        ol <- fread(file = str_c("/data/bio-tanyj/hg19/EpiMap/Enhancers/", x)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = F) %>%
        findOverlaps(., creGR)
        TMP <- rep("N", nrow(cre_all))
        TMP[unique(ol@to)] <- "Y"
        return(TMP)
    }

    timestamp()
    cat(">>>> FindOverlaps, CREs hit Enhancer_Epi.\n")
    library(furrr)
    plan(multisession, workers = NumThreds)
    HitLs <- future_map_dfc(bed_ls, HitEn, .progress = T)
    plan(sequential)
    
    # all 833 samples
    cbind(cre_all, HitLs) %>%
        saveRDS(., file = "./EpiMap_CRE_Enhancer_HitMatrix_all_biosamples_20221218.rds")

    # # only 100 representative samples.
    # cbind(cre_all, HitLs[,colnames(HitLs) %in% SelectMeta$id]) %>%
    #     saveRDS(., file = "./EpiMap_CRE_Enhancer_HitMatrix_100_biosamples_20220804.rds")
    
    # only 160 representative samples.
    cbind(cre_all, HitLs[,colnames(HitLs) %in% SelectMeta$id]) %>%
        saveRDS(., file = "./EpiMap_CRE_Enhancer_HitMatrix_160_biosamples_20221218.rds")
    
    # Chromatin states matrix -------------------------------------------
    bed_ls <- list.files(path = "/data/bio-tanyj/hg19/EpiMap/ChromHMM_states_of_833samples/",
                        pattern = "BSS.+bed.gz") %>%
        as.list()
    names(bed_ls) <- str_replace_all(bed_ls, "(^.+?)_.+", "\\1")
    
    HitChrHMM <- function(x){
        df <- fread(file = str_c("/data/bio-tanyj/hg19/EpiMap/ChromHMM_states_of_833samples/", x)) 
        ol <- df %>%
        makeGRangesFromDataFrame(., keep.extra.columns = F, seqnames.field = "V1", start.field = "V2", end.field = "V3") %>%
        findOverlaps(., creGR)
        TMP <- rep(NA, nrow(cre_all))
        TMP[ol@to] <- df$V4[ol@from]
        TMP <- factor(TMP, levels = c("TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", 
                                    "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts", 
                                    "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", 
                                    "Quies"))
        return(TMP)
    }

    timestamp()
    cat(">>>> ChromHMM states\n")
    library(furrr)
    plan(multisession, workers = NumThreds)
    HitLs <- future_map_dfc(bed_ls, HitChrHMM, .progress = T)
    plan(sequential)
    
    
    # all 833 samples
    cbind(cre_all, HitLs) %>%
        saveRDS(., file = "./EpiMap_CRE_ChromHMM_HitMatrix_all_biosamples_20221218.rds")

    # # only 100 representative samples.
    # cbind(cre_all, HitLs[,colnames(HitLs) %in% SelectMeta$id]) %>%
    #     saveRDS(., file = "./EpiMap_CRE_ChromHMM_HitMatrix_100_biosamples_20220804_1.rds")

    # only 160 representative samples.
    cbind(cre_all, HitLs[,colnames(HitLs) %in% SelectMeta$id]) %>%
        saveRDS(., file = "./EpiMap_CRE_ChromHMM_HitMatrix_160_biosamples_20221218.rds")
