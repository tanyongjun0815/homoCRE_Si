#!/bin/bash
#BSUB -J EpiEn                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 2                  ### ask for number of cores (default: 1)
#BSUB -W 10:20               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node


## resized TEs and CRE peaks, 20210510 
Rscript /scratch/2022-02-03/bioTanyj/UMI_silencer/scripts/R3_Transition/R_3.1.0_Epigenetic_enhancers.R

