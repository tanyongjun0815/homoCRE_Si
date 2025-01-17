#!/bin/bash
#BSUB -J GAT_eQTl                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 5                  ### ask for number of cores (default: 1)
#BSUB -W 48:20               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node


Rscript402 /scratch/2023-04-28/bioTanyj/UMI_silencer/scripts/R7_GWAS_eQTL/R_7.1_CRE_GWAS_eQTL_GAT.R

