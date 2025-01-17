#!/bin/bash
#BSUB -J bed                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 1                  ### ask for number of cores (default: 1)
#BSUB -W 2:00               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

##-------------------------------------------------------------------------------------##
# Prepare bed files (merge, split based on strand) for transformation (bed2bw) and calling.
#                                                        Tan,Yongjun 20210911, 20230101
#--------------------------------------------------------------------------------------##
# Use notExtend fragments of plasmid.

# cd /scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S2_alignment

# ## merge bed file------------------------------------------
# # HCT116 
#     # cDNA
#     cat Uniq_UMI_dT_Extend_si116_C5_S1.bed Uniq_UMI_NodT_si116_C5_S1.bed > Si116_cDNA_rep1.bed
#     cat Uniq_UMI_dT_Extend_si116_C6_S2.bed Uniq_UMI_NodT_si116_C6_S2.bed > Si116_cDNA_rep2.bed
#     # plasmid
#     cat Uniq_UMI_dT_notExtend_si116_P1_S3.bed Uniq_UMI_NodT_si116_P1_S3.bed > Si116_plasmid_rep1.bed
#     cat Uniq_UMI_dT_notExtend_si116_P2_S4.bed Uniq_UMI_NodT_si116_P2_S4.bed > Si116_plasmid_rep2.bed
#     # merge
#     cat Si116_cDNA_rep1.bed Si116_cDNA_rep2.bed > Si116_cDNA_Merge.bed
#     cat Si116_plasmid_rep1.bed Si116_plasmid_rep2.bed > Si116_plasmid_Merge.bed

# # K562 (two batches data of cDNA were generated.)
#     # cDNA
#     cat Uniq_UMI_dT_Extend_si562_C9_S1.bed \
#         Uniq_UMI_dT_Extend_si562_C25_S1.bed \
#         Uniq_UMI_NodT_si562_C9_S1.bed \
#         Uniq_UMI_NodT_si562_C25_S1.bed > Si562_cDNA_rep1.bed
#     cat Uniq_UMI_dT_Extend_si562_C10_S2.bed \
#         Uniq_UMI_dT_Extend_si562_C27_S1.bed \
#         Uniq_UMI_NodT_si562_C10_S2.bed \
#         Uniq_UMI_NodT_si562_C27_S1.bed > Si562_cDNA_rep2.bed
#     # plasmid
#     cat Uniq_UMI_dT_notExtend_si562_P11_S3.bed Uniq_UMI_NodT_si562_P11_S3.bed > Si562_plasmid_rep1.bed
#     cat Uniq_UMI_dT_notExtend_si562_P12_S4.bed Uniq_UMI_NodT_si562_P12_S4.bed > Si562_plasmid_rep2.bed
#     # merge
#     cat Si562_cDNA_rep1.bed Si562_cDNA_rep2.bed > Si562_cDNA_Merge.bed
#     cat Si562_plasmid_rep1.bed Si562_plasmid_rep2.bed > Si562_plasmid_Merge.bed

# ## split into each strand -------------------------------------
# cd /scratch/2021-09-09/bioTanyj/UMI_silencer/R1_identification/S2_alignment
# for i in `ls | grep "Merge.bed"`
# do  
#     {
#         echo ">>> ${i} `date`"
#         # #count
#         # echo "${i}, `cat ${i} | wc -l`" >> Number_of_fragments_generated.txt

#         # split
#         for j in + -
#         do
#             cat ${i} | grep ${j} > ${i%Merge.bed}Strand${j}.bed
#         done
#     }&
# done

## each chromosome ----------------------------------------------
cd /scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S3_calling
for i in `ls | grep -E "Strand|Merge|rep" | grep -v chr`
do
    echo ">>>${i}"
    # for j in {1..22} X Y
    for j in Y M
    do
        {
        echo "      chr${j}"
        cat ${i} | \
            awk -v CHR=chr${j} '{if($1==CHR){print $0}}' | sort -k1,1 -k2,2n -k3,3n \
            > ${i%.bed}_chr${j}.bed
        }&
    done
done

# # each replictes
# for i in `ls | grep -E "rep" | grep -v chr`
# do
#     echo ">>>${i}"
#     for j in {1..22} X Y
#     do
#         {
#         echo "      chr${j}"
#         cat ${i} | \
#             awk -v CHR=chr${j} '{if($1==CHR){print $0}}' | sort -k1,1 -k2,2n -k3,3n \
#             > ${i%.bed}_chr${j}.bed
#         }&
#     done
# done


# cd /scratch/2022-02-03/bioTanyj/UMI_silencer/R1_identification/S3_calling
# for i in `ls | grep -E ".+Merge.+" | grep -v chr`
# do
#     echo ">>>${i}"
#     for j in {1..22} X Y
#     do
#         {
#         echo "      chr${j}"
#         cat ${i} | \
#             awk -v CHR=chr${j} '{if($1==CHR){print $0}}' | sort -k1,1 -k2,2n -k3,3n \
#             > ${i%.bed}_chr${j}.bed
#         }&
#     done
# done