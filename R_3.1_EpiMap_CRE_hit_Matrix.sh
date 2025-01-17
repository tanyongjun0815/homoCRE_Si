#!/bin/bash
#BSUB -J EpiMap                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 6                  ### ask for number of cores (default: 1)
#BSUB -W 3:20               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

Rscript402 /scratch/2022-12-14/bioTanyj/UMI_silencer/scripts/R3_Transition/R_3.1_EpiMap_CRE_hit_Matrix.R