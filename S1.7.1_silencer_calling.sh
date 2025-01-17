#!/bin/bash
#BSUB -J SiCalling1                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 2                  ### ask for number of cores (default: 1)
#BSUB -W 40:20             ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

##-------------------------------------------------------------------------------------##
# Mege replicates, split .bed file of each chromosome; then identification silencers.
#                                                         use smp queue in calling
#                                                        Tan,Yongjun 20210521/20210912/202212
#--------------------------------------------------------------------------------------##

cd /scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S3_calling
referenceChrom="/data/bio-tanyj/hg19/homo_sapiens_normal.chrom.sizes"


###---------------------------------- Calling -------------------------------------
## Silencers 20210912 (each replicate/strand, merged data)
for i in `ls | grep -v chr | grep cDNA`
do
    {
    # for j in {22..1} X
    for j in Y M
    do
          plasmid_name=${i//cDNA/plasmid}
          cDNA_name=${i}

          tempfile=${plasmid_name%.bed}_chr${j}.Silencer_finding
          resultfile=${tempfile//_plasmid_/_}
          resultfile=Silencers_${resultfile%.Silencer_finding}_UMI_rmdup_padjBH.tsv
          # echo "${cDNA_name%.bed}_chr${j}.bed :: ${plasmid_name%.bed}_chr${j}.bed;  ${tempfile};  ${resultfile}"

          if [ -f "$tempfile" ]
          then
            echo "      $tempfile found."
          else
            echo ">>>>> $tempfile not found, start identify silencers.........."
            touch ${tempfile}
            referenceChrom="/data/bio-tanyj/hg19/chromosome_size/chr${j}.size"

            # run identification of silencer on this chromosome.
            Rscript402 /scratch/2023-01-01/bioTanyj/silencer_old/scripts_202011/insulatorFindUseP_adjust_modified.r \
                          control:${plasmid_name%.bed}_chr${j}.bed \
                          sample:${cDNA_name%.bed}_chr${j}.bed \
                          chrom:${referenceChrom} \
                          output:${resultfile}
          fi
    done
    }&
    # wait
done



# #### 202108
# ## calling (bed files generated using UMI_rmdup.bed, each strand)------------------
# for	i in - +
# do
#     {
#         for j in Y X {1..22}
#         do
#             tempfile=Silencer_Strand_chr${j}_si116_Strand${i}_UMI_rmdup.finding
#             if [ -f "$tempfile" ]
#             then
#               echo "      $tempfile found."
#             else
#               echo ">>>>> $tempfile not found, start identify silencers.........."
#               touch ${tempfile}
#               referenceChrom="/data/bio-tanyj/hg19/chromosome_size/chr${j}.size"

#               # run identification of silencer on this chromosome.
#               Rscript /scratch/2021-08-23/bioTanyj/silencer_old/scripts_202011/insulatorFindUseP_adjust_modified.r \
#                             control:chr${j}_Uniq_UMI_si116_PM_Strand${i}.bed \
#                             sample:chr${j}_Uniq_UMI_si116_CM_Strand${i}.bed \
#                             chrom:${referenceChrom} \
#                             output:Silencer_Strand${i}_chr${j}_UMI_rmdup_padjBH.tsv
#             fi
#         done
#     }&
# done



# ## calling (bed files generated using ..._OldMethod_rmdup.bed)------------------
# for	i in 1 2 M
# do
#     {
#         for j in Y X {1..22}
#         do
#           tempfile=Silencer_chr${j}_si116_sample${i}_OldMethod_rmdup.finding
#           if [ -f "$tempfile" ]
#           then
#             echo "      $tempfile found."
#           else
#             echo ">>>>> $tempfile not found, start identify silencers.........."
#             touch ${tempfile}
#             referenceChrom="/data/bio-tanyj/hg19/chromosome_size/chr${j}.size"

#             # run identification of silencer on this chromosome.
#             Rscript /scratch/2021-08-16/bioTanyj/silencer_old/scripts_202011/insulatorFindUseP_adjust_modified.r \
#                           control:chr${j}_Uniq_Old_strategy_si116_P${i}.bed \
#                           sample:chr${j}_Uniq_Old_strategy_si116_C${i}.bed \
#                           chrom:${referenceChrom} \
#                           output:Silencer_chr${j}_si116_sample${i}_OldMethod_rmdup_padjBH.tsv
#           fi
#         done
#     }&
# done

# ## calling (bed files generated using ..._UMI_rmdup.bed)-----------------------
# for	i in 1 2 M
# do
#     {
#         for j in Y X {1..22}
#         do
#           tempfile=Silencer_chr${j}_si116_sample${i}_UMI_rmdup.finding
#           if [ -f "$tempfile" ]
#           then
#             echo "      $tempfile found."
#           else
#             echo ">>>>> $tempfile not found, start identify silencers.........."
#             touch ${tempfile}
#             referenceChrom="/data/bio-tanyj/hg19/chromosome_size/chr${j}.size"

#             # run identification of silencer on this chromosome.
#             Rscript /scratch/2021-08-16/bioTanyj/silencer_old/scripts_202011/insulatorFindUseP_adjust_modified.r \
#                           control:chr${j}_Uniq_UMI_si116_P${i}.bed \
#                           sample:chr${j}_Uniq_UMI_si116_C${i}.bed \
#                           chrom:${referenceChrom} \
#                           output:Silencer_chr${j}_si116_sample${i}_UMI_rmdup_padjBH.tsv
#           fi
#         done
#     }&
# done

