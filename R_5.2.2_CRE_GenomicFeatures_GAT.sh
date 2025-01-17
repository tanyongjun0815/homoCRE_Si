#!/bin/bash
#BSUB -J GFeature                  ### set the job Name
#BSUB -q short               ### specify queue, medium/short/debug/ser
#BSUB -n 40                  ### ask for number of cores (default: 1)
#BSUB -W 20:20              ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node


## resized TEs and CRE peaks, 20210510 
Rscript /scratch/2022-12-20/bioTanyj/UMI_silencer/scripts/R5_Distribution/R_5.2.2_CRE_GenomicFeatures_GAT.R & \
    Rscript /scratch/2022-12-20/bioTanyj/UMI_silencer/scripts/R5_Distribution/R_5.2.2_CRE_GenomicFeatures_GAT.R & \
    Rscript /scratch/2022-12-20/bioTanyj/UMI_silencer/scripts/R5_Distribution/R_5.2.2_CRE_GenomicFeatures_GAT.R

