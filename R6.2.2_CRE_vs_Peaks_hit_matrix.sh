#!/bin/bash
#BSUB -J TFHitM                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 10                  ### ask for number of cores (default: 1)
#BSUB -W 10:00               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

Rscript /scratch/2022-12-25/bioTanyj/UMI_silencer/scripts/R6_sequence_feature/R6.2.2_CRE_vs_Peaks_hit_matrix.R