#######################################################################
# Calculate interactions between CREs and flanking regions (<3MB).
#     CREs located to other kind of CREs less than 10kb were discarded.
#                                      Tanyongjun, 20201226, 20220312
#######################################################################

  setwd("/scratch/2022-03-10/bioTanyj/UMI_silencer/R5_Distribution/Interaction_frequency")
  source("/work/bio-tanyj/soft/TanYongjun_code.R")
  library(furrr)
  plan(multicore, workers = 38)
  options(scipen = 200)
  
  cis_all <- readRDS("cre_no_heter_CREcluster_for_strawR.rds") %>%
      dplyr::filter(seqnames %in% paste("chr", c(1:22, "X"), sep = ""),
                    !is.na(ID))
  WidthFlank <- 3e6
  Binsize <- 10000
  SubsetCRE <- F # Not needed when future_map() was used! just calculate interactions of xxx CREs and corresponding flanking regions to accelerate calculation.
  NumSubsetCRE <- 40000 # Number of CREs more than 40000 will need more than 24 hours.
  # table(cis_all$ID)

# subset interactions with straw-r for each CRE. (more convenience than straw-C)
  for (i in unique(cis_all$ID)) {
    if (!file.exists(paste(i, "_straw.calculating", sep = ""))) {
      write.table(NA, file = paste(i, "_straw.calculating", sep = ""))
      cat(">>>>Calculating ", i, format(Sys.time()), "\n")

      # prepare
      cis_tmp <- cis_all %>%
        dplyr::filter(ID == i) %>% 
        dplyr::select(1:3) %>% 
        makeGRangesFromDataFrame() %>%
        GenomicRanges::resize(width = Binsize, fix = "center") %>%
        as.data.frame() %>%
        dplyr::mutate(seqnames = str_replace_all(seqnames, "chr", "")) %>%
        dplyr::rename(chromosome = seqnames) %>%
        dplyr::select(1:3) %>%
        mutate(leftBorder = start - WidthFlank, 
               rightBorder = end + WidthFlank) %>%
        mutate(leftBorder = ifelse(leftBorder > 0, leftBorder, 0),
               id = str_c(chromosome, start, end, sep = "_")) %>%
        split(., .$id)
      
      hicFile <- ifelse(str_detect(i, "HCT116"), 
                        "/data/bio-tanyj/hic/GSE104333_HCT116_Rao-2017-untreated_combined_30.hic",
                        "/data/bio-tanyj/hic/GSE63525_K562_combined_30.hic")
      
      cat("     total number of CREs:", length(cis_tmp), "\n")

      # subset
      if (SubsetCRE & length(cis_tmp) > NumSubsetCRE) {
        cat("    Number of CREs in ", i, " was more than ", NumSubsetCRE, 
            ", just calculate interaction of  ", NumSubsetCRE, ".\n")
        set.seed(123456)
        cis_tmp <- cis_tmp[, sample(seq_len(ncol(cis_tmp)), size = NumSubsetCRE, replace = F)]
      }
      
      # calculate
      resultsTMP <- future_map_dfr(cis_tmp, function(x){
        strawr::straw("observed",
                      "KR",
                      hicFile,
                      paste(x[1], ":", x[2], ":",
                            x[3], sep = ""),
                      paste(x[1], ":", x[4],
                            ":", x[5], sep = ""),
                      "BP",
                      Binsize) %>%
          return()
      }, .id = "id", .progress = T)
      
      # merge results and calculate distance
      resultsTMP <- resultsTMP %>%
        fwrite(., file = paste(i, "_straw_extracted_clean.tab", sep = ""))
        
      cat("------", i, " finished", format(Sys.time()), "\n")
    } # if
  } # i

  plan(sequential)