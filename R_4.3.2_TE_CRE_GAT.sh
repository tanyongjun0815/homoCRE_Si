#!/bin/bash
#BSUB -J GAT_TEs                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 2                  ### ask for number of cores (default: 1)
#BSUB -W 100:20               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node


## resized TEs and CRE peaks, 20210510 
Rscript402 /scratch/2023-05-29/bioTanyj/UMI_silencer/scripts/R4_Conservation/R_4.3.2_TE_CRE_GAT.R 

