#!/bin/bash
#BSUB -J QC                  ### set the job Name
#BSUB -q ser               ### specify queue, medium/short/debug/ser
#BSUB -n 16                  ### ask for number of cores (default: 1)
#BSUB -W 50:3             ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

# ## all raw reads -----------------
# cd /scratch/2021-05-13/bioTanyj/UMI_STARR_HCT116_Si/PM-XS01KF2020070275-28_HCT116_silencer_UMI_202105/PM-XS01KF2020070275-28_basecalling

# /work/bio-tanyj/soft/fastqc_v0.11.9/FastQC/fastqc -t 15 *fastq.gz -o ./QC/

# cd ./QC
# /work/bio-tanyj/miniconda3/bin/multiqc ./

# ## dT and nodT reads (20210817)--------------------
# cd /scratch/2021-08-16/bioTanyj/UMI_silencer/UMI_STARR_HCT116_Si/PM-XS01KF2020070275-28_HCT116_silencer_UMI_202105/alignment/fastq_trimmed

# /work/bio-tanyj/soft/fastqc_v0.11.9/FastQC/fastqc -t 16 *fastq.gz -o ./

# /work/bio-tanyj/miniconda3/bin/multiqc ./

## dT and nodT reads (20210817)--------------------
cd /scratch/2021-08-16/bioTanyj/UMI_silencer/UMI_STARR_HCT116_Si/PM-XS01KF2020070275-28_HCT116_silencer_UMI_202105/alignment/fastq_trimmed_0bp_mismatched

/work/bio-tanyj/soft/fastqc_v0.11.9/FastQC/fastqc -t 16 *fastq.gz -o ./

/work/bio-tanyj/miniconda3/bin/multiqc ./

