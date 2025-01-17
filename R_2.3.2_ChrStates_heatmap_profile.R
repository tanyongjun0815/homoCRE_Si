###########################################################################################
#    calculate reference-point based score of chromatin 
#            states (.bw file) of CRE in .bed files.
#                                           Tanyongjun 20210314/20210321/20211023/20211210
###########################################################################################

setwd("/scratch/2021-12-07/bioTanyj/UMI_silencer/R2_ChrState/enrichment_profile")
source("/work/bio-tanyj/soft/TanYongjun_code.R")
library(stringr)
NumThreads <- 40

system("source activate base")

#------------------------------20211023------------------------------------
OnlyChrState <- T # TRUE: only analysis chromatin states realted (DHS, H2AFZ, H3K4me.....)

cre_type_ls <- list.files(path = "./bed_files/", pattern = "*.bed") %>%
    str_replace_all(., ".+?_(.+?)_.+", "\\1") %>%
    str_replace_all(., "\\d+", "") %>%
    str_replace_all(., "\\.|\\+|-", "") %>%
    str_replace_all(., c("CellSpeciCellSpeci" = "CellSpeci",
                         "CellSpeciCommon" = "CellSpeci")) %>%
    str_replace_all(., "noDHS", "DHS") %>%
    unique()
for(cre_type in cre_type_ls){ # ALL, Kmenas, ActivityCluster, Trans, Cons.
    bed_ls <- list.files(path = "./bed_files/", pattern = str_c(".+", cre_type, ".+", sep = ""))
    cre_id_ls <- str_replace_all(bed_ls, "(.+?)_(.+?)_.+", "\\1") %>%
        unique()
    
    # create dir for each cre_type
    system(command = str_c("mkdir ", cre_type), ignore.stderr = T)
    setwd(str_c("./", cre_type))

    # for ALL (each kind of CRE as a whole), plot all kind CRE in one figure.
    if (cre_type == "All") {
       cre_id_ls = "All"
    }

    for(cre_id in cre_id_ls){ # Silencer, enhancer, promoter, Random, and All.
        name_cre_bed <- bed_ls[str_detect(bed_ls, cre_id)]
        
        for(cell_line in c("HCT116", "K562")){ # cell line
            name_cre_bed_cl <- name_cre_bed[str_detect(name_cre_bed, cell_line)]
            cat(">>>CRE type:", cre_type, 
                "  \nCRE id:", cre_id, 
                "  \ncell line:", cell_line,
                "  \n", str_c(name_cre_bed_cl, collapse = "\n"), 
                "\n\n")
            
            # list of bw file
            switch(cell_line, 
                   "HCT116" = bwFile_path <- "/data/bio-tanyj/hg19/encode/hct116/",
                   "K562" = bwFile_path <- "/data/bio-tanyj/hg19/encode/k562/")
            bwFile_ls <- list.files(path = bwFile_path, pattern = str_c(cell_line, "_.+.bigWig", sep = ""))

            if (OnlyChrState) { # only chromatin state-realted datasets?
               bwFile_ls <- bwFile_ls[str_detect(bwFile_ls, "H2AZ|DHS|H3|H4|DNAmethy")]
            }
            # 
            for(bwFile in bwFile_ls){# each bw file.
                cat(">>> ", bwFile, "\n")
                id <- str_c(cell_line, cre_type, cre_id, 
                            str_replace_all(bwFile, "(.+)_hg19.+", "\\1"), sep = "_")
                id_lock <- str_c(id, ".running")
                
                if (!file.exists(id_lock)) {
                    write.table(NA, id_lock)

                    # calculate score
                    system(paste("/work/bio-tanyj/miniconda3/bin/computeMatrix reference-point ",
                                " --referencePoint center ",
                                " -S ", bwFile_path, bwFile, " ",
                                " -R ", str_c(str_c("../bed_files/", name_cre_bed_cl, sep = ""), collapse = " "), " ",
                                " -o ", id, "_computematrix.gz ",
                #                 "--outFileNameMatrix ", id, "_computematrix.matrix ",
                #                 "--outFileSortedRegions ", id, "_computematrix.bed ",
                                " -b 5000 -a 5000 --binSize 20 -p ", NumThreads,
                                " --smartLabels",
                                sep = ""))

                    # plot profile
                    system(paste("/work/bio-tanyj/miniconda3/bin/plotProfile ",
                                "-m ", id, "_computematrix.gz ",
                                "-out Profile_", id, ".png ",
                                "--outFileNameData ", id, "_plotProfile.tab ",
                                #  "--plotTitle ", id, " ",
                                "--plotWidth 16 --plotHeight 12 ",
                                sep = ""), wait = F)

                    # plot heatmap
                    system(paste("/work/bio-tanyj/miniconda3/bin/plotHeatmap ",
                                "-m ", id, "_computematrix.gz ",
                                "-out Heatmap_", id, ".png ",
                                "--colorMap gist_earth ",
                                #  "--outFileNameMatrix PlotHeatmap_matrix_", id, ".gz ",
                                "--dpi 300 ",
                                "--heatmapHeight 38 --heatmapWidth 16",
                                sep = ""), wait = F)

                } # if
            } # bw file
        } # cell_line
    } # cre_id

    # return uplevel dir
    setwd("../")
}# cre_type

##------------------- check -------------------------------------------
for(i in list.files(pattern = ".+running", recursive = T)){
    profile_name <- str_replace_all(i, ".running", ".png") %>%
        str_replace_all(., "/", "/Profile_")
    if (file.exists(str_c("./", profile_name, sep = ""))) {
        # cat(" ")
    }else{
        cat(">>> ", profile_name, "not exist!\n")
    }
}



# ##--------------------------------------- A.t code -------------------
# ## compute matrix and plot heatmap
# for (i in list.files(path = "../bw_file_chromatin_states", pattern = ".+.bw$")){
#     # i <- "DNAmethy.bw"
#     id <- str_replace_all(i, ".bw", "_10kbWidth_10bpBin")
#     # compute reference-point matrix
#     if(!file.exists(paste(id, "_computematrix.calculating", sep = ""))){
#         timestamp()
#         cat(">>> ", id, "\n")
#         write.table(id, file = paste(id, "_computematrix.calculating", sep = ""))
#         system(paste("/work/bio-tanyj/miniconda3/bin/computeMatrix reference-point ",
#                     "--referencePoint center ",
#                     "-S ../bw_file_chromatin_states/",
#                     i, " ",
#                     "-R ",
#                     str_c(unique(cis_all$ID), collapse = " "), " ",
#                     "-o ", id, "_computematrix.gz ",
#     #                 "--outFileNameMatrix ", id, "_computematrix.matrix ",
#     #                 "--outFileSortedRegions ", id, "_computematrix.bed ",
#                     "-b 5000 -a 5000 --binSize 10 -p ", NumThreads,
#                     sep = ""))

#         # plot profile
#         system(paste("/work/bio-tanyj/miniconda3/bin/plotProfile ",
#                     "-m ", id, "_computematrix.gz ",
#                     "-out Profile_", id, ".png ",
#                     "--outFileNameData ", id, "_plotProfile.tab ",
#                     #  "--plotTitle ", id, " ",
#                     "--plotWidth 16 --plotHeight 12 ",
#                     sep = ""))

#         # plot heatmap
#         system(paste("/work/bio-tanyj/miniconda3/bin/plotHeatmap ",
#                     "-m ", id, "_computematrix.gz ",
#                     "-out Heatmap_", id, "-1.png ",
#                     "--colorMap gist_earth ",
#                     #  "--outFileNameMatrix PlotHeatmap_matrix_", id, ".gz ",
#                     "--dpi 300 ",
#                     "--heatmapHeight 38 --heatmapWidth 16",
#                     sep = ""))
#     } # if
# } # for

# ## plot
# for (i in list.files(path = "../bw_file_chromatin_states", pattern = ".+.bw$")){
#     # i <- "DNAmethy.bw"
#     id <- str_replace_all(i, ".bw", "_10kbWidth_10bpBin")
#     # compute reference-point matrix
#     if(!file.exists(paste(id, "_computematrix.ploting", sep = ""))){
#         timestamp()
#         cat(">>> ", id, "\n")
#         write.table(id, file = paste(id, "_computematrix.ploting", sep = ""))

#         # plot profile
#         # system(paste("/work/bio-tanyj/miniconda3/bin/plotProfile ",
#         #             "-m ", id, "_computematrix.gz ",
#         #             "-out Profile_", id, ".png ",
#         #             "--outFileNameData ", id, "_plotProfile.tab ",
#         #             #  "--plotTitle ", id, " ",
#         #             "--plotWidth 16 --plotHeight 12 ",
#         #             sep = ""))

#         # plot heatmap
#         system(paste("/work/bio-tanyj/miniconda3/bin/plotHeatmap ",
#                     "-m ", id, "_computematrix.gz ",
#                     "-out Heatmap_", id, "_x.png ",
#                     "--colorMap gist_earth ",
#                     #  "--outFileNameMatrix PlotHeatmap_matrix_", id, ".gz ",
#                     "--dpi 300  --perGroup ",
#                     "--heatmapHeight 16 --heatmapWidth 4",
#                     sep = ""))
#     } # if
# } # for
