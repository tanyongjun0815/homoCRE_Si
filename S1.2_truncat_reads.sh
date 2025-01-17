#!/bin/bash
#BSUB -J truncate                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 5                  ### ask for number of cores (default: 1)
#BSUB -W 35:3             ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

##################################################################################
#           Truncate forward and reverse reads.


## test ------------------------------------------
# truncate with seqkit (test data)
# cd /scratch/2021-05-13/bioTanyj/UMI_STARR_HCT116_Si/PM-XS01KF2020070275-28_HCT116_silencer_UMI_202105/PM-XS01KF2020070275-28_basecalling/test
#     echo ">>> truncate"
#     for j in 1 3
#     do
        
#         {
#             /work/bio-tanyj/miniconda3/bin/seqtk trimfq -e 55 -b 55 Test_si116_P1_S3_L002_R${j}_001.fastq.gz | gzip > Test_si116_P1_S3_L002_R${j}_trimmed.fastq.gz
#         }&
#     done
#     wait

## real ------------------------------------------

# # K562 batch1 truncate with seqtk (36bp) 20210910 ----------------------
# cd /scratch/2023-01-01/bioTanyj/UMI_silencer/UMI_STARR_K562_Si/basecalling/batch1
# CleanData="/scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S1_truncate"
# Trim5end="25" # cut xx bp from 5'end
# Trim3end="89" # cut xx bp from 3'end

# for i in `ls | grep -E "si562.+001.fastq.gz"  | grep -v "R2"`
# do
#     {
#         /work/bio-tanyj/miniconda3/bin/seqtk trimfq -b ${Trim5end} -e ${Trim3end} ${i} | gzip > ${CleanData}/${i%_001.fastq.gz}_trim36bp.fastq.gz
#     }&
# done
# wait

# # K562 batch2 truncate with seqtk (36bp) 20210910 ----------------------
# cd /scratch/2023-01-01/bioTanyj/UMI_silencer/UMI_STARR_K562_Si/basecalling/batch2
# CleanData="/scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S1_truncate"

# for i in `ls | grep -E "si562.+001.fastq.gz"  | grep -v "R2"`
# do
#     {
#         /work/bio-tanyj/miniconda3/bin/seqtk trimfq -b ${Trim5end} -e ${Trim3end} ${i} | gzip > ${CleanData}/${i%_001.fastq.gz}_trim36bp.fastq.gz
#     }&
# done
# wait

# HCT116 truncate with seqtk (36bp) 20210910 ----------------------
cd /scratch/2023-01-01/bioTanyj/UMI_silencer/UMI_STARR_HCT116_Si/PM-XS01KF2020070275-28_HCT116_silencer_UMI_202105/PM-XS01KF2020070275-28_basecalling
CleanData="/scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S1_truncate"

for i in `ls | grep -E "si116.+001.fastq.gz"  | grep -v "R2"`
do
    {
        /work/bio-tanyj/miniconda3/bin/seqtk trimfq -b ${Trim5end} -e ${Trim3end} ${i} | gzip > ${CleanData}/${i%_001.fastq.gz}_trim36bp.fastq.gz
    }&
done
wait