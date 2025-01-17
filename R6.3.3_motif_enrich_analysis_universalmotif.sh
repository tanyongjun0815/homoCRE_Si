#!/bin/bash
#BSUB -J UniEnri                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 15                  ### ask for number of cores (default: 1)
#BSUB -W 48:20               ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

source activate UMI
/work/bio-tanyj/miniconda3/envs/UMI/bin/Rscript /scratch/2023-04-18/bioTanyj/UMI_silencer/scripts/R6_sequence_feature/R6.3.3_motif_enrich_analysis_universalmotif.R