#!/bin/bash
#BSUB -J InterAct                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 20                  ### ask for number of cores (default: 1)
#BSUB -W 24:00               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

cd /scratch/2022-03-10/bioTanyj/UMI_silencer/R5_Distribution/Interaction_frequency

Rscript /scratch/2022-03-10/bioTanyj/UMI_silencer/scripts/R5_Distribution/R_5.4.3_extract_interactions_from_hic_with_strawR_clean.R

