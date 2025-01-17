#!/bin/bash
#BSUB -J TFsHis                  ### set the job Name
#BSUB -q short               ### specify queue, medium/short/debug/ser
#BSUB -n 40                 ### ask for number of cores (default: 1)
#BSUB -W 100:5               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

# when ask 40 cores, only 50%-60% of CPU was used. Thus, 20 cores was a better choce.

Rscript /scratch/2022-09-21/bioTanyj/UMI_silencer/scripts/R2_ChrState/R_2.2.1_CRE_TF_histone_GAT.R
