##################################################################################
#  subset population index in CRE; evenly subset from whole genome for plot on PC.
#                         Tanyongjun 20210418 modified 20210416
##################################################################################

rm(list = ls())
setwd("/scratch/2022-03-10/bioTanyj/Hapmap_1000Genomes/PopIndex")
source("/work/bio-tanyj/soft/TanYongjun_code.R")

#### merge -------------------------------------------------
# load pi (ref genome: GRCh38)
    pi_all <- NULL
    for(i in list.files(path = "./", pattern = ".+windowed.pi$")){
      cat(">   ", i, "\n")
        pi_all <- fread(i) %>%
          dplyr::rename(seqnames = CHROM, 
                        start = BIN_START, end = BIN_END) %>%
          rbind(., pi_all)
    }
    cat("\n>>>>>>> load pi completed!\n")

    pi_all %>%
      slice_sample(., n = 1000000) %>%
      fwrite(., file = "./Pi_all_sampled_1000000_hg38_regions.csv.gz")

# hg38 to hg19
    cat("\n>>>>>>> liftOverl hg38 to hg19!\n")
    library(liftOver)
    ch1 = import.chain("/scratch/2022-03-10/bioTanyj/UMI_silencer/data_file/liftOver_chainfile/hg38ToHg19.over.chain")
    pi_gr19 <- pi_all %>%
      # head(., 200000) %>%
      dplyr::mutate(queryHits = 1:nrow(.)) %>%
      makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
      GenomicRanges::resize(width = 50, fix = "center") %>%
      liftOver(., ch1) %>% 
      unlist()
    pi_dfhg19 <- as.data.frame(pi_gr19) %>%
      dplyr::select(queryHits, PI)


#### CRE score -------------------

    cre_all <- read_rds("/scratch/2022-03-10/bioTanyj/UMI_silencer/data_file/SiMerge_EnhancerENCODE_SiReported_Pro_Random_KmeansCluster7_TransCluster4_ConsCluster_20211025.rds") %>%
      # slice_sample(., n = 100000) %>%
      dplyr::mutate(subjectHits = 1:nrow(.))

    cre_allGR <- makeGRangesFromDataFrame(cre_all, keep.extra.columns = F) %>%
      GenomicRanges::resize(., width = 5, fix = "center")

    ol <- findOverlaps(pi_gr19, cre_allGR) %>%
      as.data.frame() %>%
      left_join(cre_all, .) %>%
      left_join(., pi_dfhg19) %>%
      dplyr::select(-subjectHits, -queryHits)

    ol_mean <- ol %>%
      dplyr::group_by(seqnames, start, end, CRE, cell_line, sample, ID) %>%
      dplyr::summarise(PI = mean(PI, na.rm = T), N_variants = )

    fwrite(ol_mean, "PI_CRE_all_mean_20220308.tsv.gz", sep = "\t")

    ol_mean %>%
      group_by(ID) %>%
      dplyr::summarise(NoMissRate = mean(!is.na(PI)))
