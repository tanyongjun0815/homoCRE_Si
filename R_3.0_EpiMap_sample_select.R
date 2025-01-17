setwd("/data/bio-tanyj/hg19/EpiMap")

## load DHS matrix and subset.
if(F){
    library(rhdf5)
    h5ls("DNase-seq_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5")
    DHS_matrix <- h5read("DNase-seq_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5", "matrix")
    # DHS_matrix = H5Fopen(name = "DNase-seq_all_bin_dense_on_mixed_impobs_r25_e100_allchr_merged.hdf5") # error reported!

    cat(dim(DHS_matrix))

    set.seed(12345)
    DHS_matrix_sub <- DHS_matrix[, sample(1:ncol(DHS_matrix), size = 10000, replace = F)]
    saveRDS(DHS_matrix_sub, file = "DHS_matrix_10000columns.rds")

    DHS_matrix_sub <- DHS_matrix[, sample(1:ncol(DHS_matrix), size = 20000, replace = F)]
    saveRDS(DHS_matrix_sub, file = "DHS_matrix_20000columns.rds")
}

# Hierarchical clustering samples.
    M_rowname <- read.table("./mark_matrix_names.txt")
    exp_m <- readRDS("DHS_matrix_20000columns.rds")
    row.names(exp_m) <- M_rowname$V1

    d=dist(t(exp_m))
    hc=hclust(d,method="complete")

    plot(hc,labels = FALSE, hang= -1)
    rect.hclust(hc, h = 40, border = "red")

    clu.k5=cutree(hc,k=100)
    clu.h80=cutree(hcl,h=80)
    NA  