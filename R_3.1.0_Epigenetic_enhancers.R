################################################################################
# Two steps related to identify ActiveEnhancer and PoisedEnhancer.
#   1. merge epigenetic datasets of HCT116 and K562, generate new summary file.
#   2. identify ActiveEnhancer and PoisedEnhancers.
#                               TanYongjun, 20220208
#################################################################################


# ######################################## Part 1 ######################################
# setwd("/data/bio-tanyj/hg19/histone_other_tissue")
# source("/work/bio-tanyj/soft/TanYongjun_code.R")
# ### summary of download bed files.
#   df <- tibble()
#   fileLS <- list.files(pattern = ".+hg19.+bed.gz$")
#   fileLS <- fileLS[str_detect(fileLS, "HCT116|K562")]
#   for(i in fileLS) {
#     cat(i, "\n")
#     tmp <- fread(i) %>%
#       dplyr::select(1:3) %>%
#       GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1",
#                                               start.field = "V2",
#                                               end.field = "V3")
#     df <- rbind(df,
#                 c(i, length(tmp), sum(width(tmp)) / 1e6,
#                 quantile(width(tmp), probs = seq(0, 1, 0.1))))
#   }
#   colnames(df) <- c("ID", "Num_regions", "TotalLength(MB)",
#                 paste(0:10 * 10, "%", sep = ""))

#   fwrite(df, file = "Summary_of_download_bed_files_20220208.csv")


# ### summary of each cell x feature (calculate the total length, generate merged file)

# allFile <- data.frame(ID = list.files(pattern = ".+hg19.+bed.gz$"), 
#                       cell_line = NA, 
#                       features = NA)
# allFile <- allFile %>%
#   dplyr::mutate(cell_line = str_replace_all(ID, 
#                                             "(.+)_(DHS|H2AFZ|H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K79me2|H3K9ac|H3K9me2|H3K9me3|H4K20me1)_hg19.+",
#                                             "\\1"),
#                 features = str_replace_all(ID, 
#                                             "(.+)_(DHS|H2AFZ|H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K79me2|H3K9ac|H3K9me2|H3K9me3|H4K20me1)_hg19.+",
#                                             "\\2")) %>%
#   dplyr::filter(str_detect(cell_line, "HCT116|K562")) # only re-generate merged files of HCT116 and K562.

# df <- NULL
# for(i in unique(allFile$cell_line)) {
#   for(j in unique(allFile$features)) {
#     FileList <- allFile$ID[allFile$cell_line == i &
#                              allFile$features == j &
#                              str_detect(allFile$ID, "merged", negate = T)]
#     if(length(FileList) > 0) {
#       cat(i, j, " Number of bed files:", length(FileList), "Total covered region length:")
#       tmp <- NULL
#       for(k in FileList) {
#         tmp <- fread(k) %>%
#                   dplyr::select(1:3) %>%
#                   rbind(tmp, .)
#       }
#       names(tmp)[1:3] <- c("chromosome", "start", "end")
#       tmp <- tmp %>%
#         dplyr::filter(chromosome != "chrY") %>% 
#         GenomicRanges::makeGRangesFromDataFrame() %>%
#         GenomicRanges::reduce()
#       cat(sum(width(tmp)) / 1e6, "\n")
#       df <- rbind(df, c(i, j, length(FileList), sum(width(tmp)) / 1e6,
#                         quantile(width(tmp), probs = seq(0, 1, 0.1))))
#       fwrite(as.data.frame(tmp), file = paste(i, j, "merged.bed.gz", sep = "_"))
#     }
#   }
# }

# colnames(df) <- c("cell_line", "features", "NumberFiles", "TotalLength(MB)",
#                 paste(0:10 * 10, "%", sep = ""))

# fwrite(df, file = "Summary_of_features_of_each_cell_line_20220208.csv")


# ######################################## Part 2 ######################################
### Identify activeenhancer and Poised enhancers in other tissues and cell lines.
  setwd("/data/bio-tanyj/hg19/histone_other_tissue")
  source("/work/bio-tanyj/soft/TanYongjun_code.R")

  allFile <- data.frame(ID = list.files(pattern = ".+merged.bed.gz$"), 
                        cell_line = NA,
                        features = NA)
  allFile <- allFile %>%
    dplyr::mutate(cell_line = str_replace_all(ID, 
                                              "(.+)_(DHS|H2AFZ|H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K79me2|H3K9ac|H3K9me2|H3K9me3|H4K20me1)_.+",
                                              "\\1"),
                  features = str_replace_all(ID, 
                                              "(.+)_(DHS|H2AFZ|H3K27ac|H3K27me3|H3K36me3|H3K4me1|H3K4me2|H3K4me3|H3K79me2|H3K9ac|H3K9me2|H3K9me3|H4K20me1)_.+",
                                              "\\2"))
  # table(allFile$cell_line, allFile$features)
  # unique(allFile$cell_line)
  # unique(allFile$features)

## identify Active/Poised enhancer.
  ReadToGR <- function(file) {
    fread(file) %>%
      makeGRangesFromDataFrame() %>%
      return(.)
  }

  enhancerALL <- list()
  allFile <- allFile %>%
    dplyr::filter(str_detect(features, "(DHS|H3K27ac|H3K27me3|H3K4me1)"))
  for(i in unique(allFile$cell_line)) {
      cat(">>> Identification enhancer of ", i, ":\n")
      cat("     Loading files:\n")
      dhsTMP <- ReadToGR(paste(i, "_DHS_merged.bed.gz", sep = ""))
      H3K27acTMP <- ReadToGR(paste(i, "_H3K27ac_merged.bed.gz", sep = ""))
      H3K4me1TMP <- ReadToGR(paste(i, "_H3K4me1_merged.bed.gz", sep = ""))

      # Active enhancers (DHS + H3K4me1 + H3K27ac)
      cat("     Active enhancer.....; ")
      activeETMP <- GenomicRanges::intersect(dhsTMP, H3K4me1TMP) %>%
        GenomicRanges::intersect(., H3K27acTMP) %>%
        # GenomicRanges::resize(width = 600, fix = "center") %>%
        as.data.frame()
      activeETMP$cell_line <- i
      activeETMP$cis_element <- "ActiveEnhancer"

      # Poised enhancers (DHS + H3K4me1, not H3K27ac)
      cat("Poised enhancer....;")
      max_dist <- 10000 # discard poised enhancer located less than max_dist from active enhancer
      tmp <- activeETMP %>%
        makeGRangesFromDataFrame() %>%
        GenomicRanges::resize(width = max_dist * 2, fix = "center")
      poisedEnhancer <- GenomicRanges::intersect(dhsTMP, H3K4me1TMP) # %>%
        # GenomicRanges::resize(width = 600, fix = "center")
      # discard posied enhancers located near active enhancers
      poisedEnhancer <- poisedEnhancer[-unique(findOverlaps(poisedEnhancer,tmp)@from)] %>%
        as.data.frame()
      poisedEnhancer$cell_line <- i
      poisedEnhancer$cis_element <- "PoisedEnhancer"

      # merge and save
      enhancerALL <- c(enhancerALL, list(rbind(activeETMP, poisedEnhancer)))
      names(enhancerALL)[length(enhancerALL)] <- i
      cat("complete; ", format(Sys.time()), "\n")
  }

# save to file.
  enhancerALL <- enhancerALL %>%
    purrr::reduce(., rbind)
  table(enhancerALL$cell_line, enhancerALL$cis_element)
  save(enhancerALL, file = "All_identified_enhancer_20220208.bin")
  enhancerALL %>%
    fwrite(., file = "All_identified_enhancer_20220208.csv.gz")