#!/bin/bash
#BSUB -J UMIfq2bed                  ### set the job Name
#BSUB -q smp               ### specify queue, medium/short/debug/ser
#BSUB -n 1                  ### ask for number of cores (default: 1)
#BSUB -W 24:20             ### set walltime limit: hh:mm
#BSUB -e %J.err              ### -o and -e mean append, -oo and -eo mean overwrite
#BSUB -o %J.out              ### Specify the output and error file. %J is the job-id
#BSUB -R "span[ptile=40]"    ### ask for 40 cores per node

##-------------------------------------------------------------------------------------##
# generate bed and bw files of UMI-starr from fastq.
#   Steps:
#      1. alignment to genome: .bam
#      2. bam 2 bed: raw.bed
#      3. Add UMI to bed: UMI.bed
# Notes: Quality trimming is not needed for bwa-mem, BWA-backtrack requires reads to be 
#        mapped in full length. Low-quality tail may greatly affect its sensitivity. 
#        Bwa-mem largely does local alignment. If a tail cannot be mapped well, 
#        it will be soft clipped. But the soft-clip may affect MAPQ.   ---Li,Heng
#
# Attention: less paired-mapped reads were found when using bwa instead of bowtie in mapping (36bp).
#            Number of properly mapped reads signicantly decreased when using 132bp in mapping (bwa).
#            used strategy: 36bp reads + bowtie 20210909
#----------------------------------------------------------------------------------------

###--------------------------------- paths and dirs -----------------------------------###
source /work/bio-tanyj/scripts/path_and_dir.sh
BowtieIndex="/data/bio-tanyj/hg19/homo_normal_with_YM"
# BowtiePath="/work/bio-tanyj/soft/bowtie-1.2.3-linux-x86_64"
# Bowtie2Path="/work/bio-tanyj/soft/bowtie2-2.4.1-linux-x86_64"
referenceFile="/data/bio-tanyj/hg19/homo_normal_with_YM.fa"
# fastqcPath="/work/bio-tanyj/soft/fastqc_v0.11.9/FastQC"
samtoolsPath="/work/bio-tanyj/soft/samtools-1.10"
bedtools2Path="/work/bio-tanyj/soft"
bwaPath="/work/bio-tanyj/miniconda3/bin"
bowtie2Path="/work/bio-tanyj/soft/bowtie2-2.4.1-linux-x86_64"

###----------------------------------set parameters ----------------------------###
resultPath="/scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S2_alignment"
CleanDataRoot="/scratch/2023-01-01/bioTanyj/UMI_silencer/R1_identification/S1_truncate"
cd ${resultPath}
nthreads=40

# bowtie parameter
    bowtiePara="-p ${nthreads} -X 2000 -v 3 -m 1 --best --strata --quiet --sam"
    #     -p(threads), -X (maximum insert size), -q (print nothing besides alighments), 
    #     -v no more than v mismatch -m just keep reads which matched no more than m position!!!
    #        (This mode only find in bowtie. Bowtie2 will random report one from all best matched position).
    #     --best reported the alignments in best-to-worse order 
    #     --strata just keep the best matched alignments (specified with --best), 
    #     --quiet print nothing besides alighments. --sam
# Samtools view parameter (for bwa pipeline, sam generate by bowtie filtered in alignment step)
    samtoolsViewPara="-bS -F 3852 -f 2 -q 40 -h"
    # -bS: : use sam as input, and bam as output.
    # -F 3852: discard reads based on the FLAG field, read unmapped/mate unmapped/
    #     not primary alignment/fail qual check/dup/supplemental align. (SAM FLAG translator)
    # -f 2: requires FLAG 2, properly alignmened
    # -q 40: discard MAPQ < 40
    #                  This parameter used in the STARRPeaker pipeline.


######################################  (bowtie 36bp)###################################

####-----------------------------Step1. alignment---------------------------------------------------
for i in `ls ${CleanDataRoot} | grep -E "trim36bp.fastq.gz"`
do
    id=${i%_L00*}
    if [ -f ${id}.calculating ]
    then
            echo ">>>>>> ${id} file exist!"   
    else
            echo "NA" > ${id}.calculating
            # i="si116_C5_S1_L002_R1_trim36bp.fastq.gz"
            echo ">>>>>> ${i}    ${id}    `date`"

        ##  (short, 40 cores)
        # # alignment
        #     echo ">>>>>> Mapping to reference genome.......`date`"
        #     ${BowtiePath}/bowtie ${bowtiePara} ${BowtieIndex} \
        #         -1 ${CleanDataRoot}/${i} \
        #         -2 ${CleanDataRoot}/${i%1_trim36bp.fastq.gz}3_trim36bp.fastq.gz \
        #         -S  | ${samtoolsPath}/samtools view -@ $nthreads -Sb - -o ${resultPath}/${id}.bam
        #     ${samtoolsPath}/samtools sort -n -@ ${nthreads} ${resultPath}/${id}.bam -o ${resultPath}/${id}_sort.bam
        #     ${samtoolsPath}/samtools fixmate -@ ${nthreads} ${resultPath}/${id}_sort.bam ${resultPath}/${id}_fixmate.bam
        #     # ${samtoolsPath}/samtools index -@ ${nthreads} ${resultPath}/${id}_fixmate.bam
        #     rm ${resultPath}/${id}_sort.bam ${nthreads} ${resultPath}/${id}.bam

        ## (smp, 1cores) 
        # bam to bed (Illumina read ID: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  <read>:<is filtered>:<control number>:<index>); 

            #              only keep <tile>:<x-pos>:<y-pos> of the read ID
            # columns in output: 1. chr, 2. start(read1), 3.end(read2), 4.trimmed readID (shared by R1,R2, and UMI), 5. strand (read1)
            echo ">>>>>> Convert to bed file.......`date`"
            ${bedtools2Path}/bedtools bamtobed -i ${resultPath}/${id}_fixmate.bam -mate1 -bedpe | \
                awk '{if($1==$4&&$1!="."&&$2!=0&&$5!=0){if($2>=$5&&$9=="-"&&$10=="+"){start=$5;end=$3;strand="-"}else if($2<$5&&$9=="+"&&$10=="-"){start=$2;end=$6;strand="+"};sub("^[^:]+:[^:]+:[^:]+:[^:]+:", "", $7);print $1, start, end, $7, strand}}' OFS='\t' \
                > ${resultPath}/${id}.bedpe # the second awd seen not needed!

        # add UMI to bed (only need one CPU, but ~ 140GB/155GB RAM!)
            # UMItab: tsv file contain two columns: seqID, UMI sequences.
            echo ">>>>>> Add UMI in bed file.......`date`"
            cat ${id}.bedpe | awk -v UMItab=${id}.UMItab -v ID=${id} -v OFS='\t' 'BEGIN{while(getline<UMItab){vec[$1]=$2}}{if($4 in vec){print $1, $2, $3, $5, $4, vec[$4] > ID"_UMI.bed"}else{print $0 > ID"_noUMI.bed"}}'
    fi
done

#####################################################################################################
#  perform fastq2bed using alignment and filtration parameters used by the STARRPEAKER pipeline.
#           less fragments were generated use this pipeline. this may caused by large genetic 
#                difference between Niu and the reference genome and the strigent filtration threshold.
#  The affect of read length and software on mapping rate was estimated in S1.4.
#                                         this pipeline was not used out study
#                                                                               202108
######################################################################################################

# ###################################### (bwa)###################################

# ####-----------------------------------------36bp---------------------------------------------------
# ####-----------------------------Step1. alignment---------------------------------------------------
# for i in `ls ${CleanDataRoot} | grep -E "R1_trim36bp.fastq.gz"`
# do
#     id=${i%_S*_L002*}
#     if [ -f ${id}.calculating ]
#     then
#         echo ">>>>>> ${id} file exist!"   
#     else
#         echo "NA" > ${id}.calculating
#         # i="si116_C5_S1_L002_R1_trim36bp.fastq.gz"
#         echo ">>>>>> ${i}    ${id}    `date`"
#         # alignment
#             echo ">>>>>> Mapping to reference genome.......`date`"
#             ${bwaPath}/bwa mem -t ${nthreads} ${referenceFile} \
#                 ${CleanDataRoot}/${i} \
#                 ${CleanDataRoot}/${i%R1_trim36bp.fastq.gz}R3_trim36bp.fastq.gz | \
#                 # ${samtoolsPath}/samtools sort -n -@ ${nthreads} | \
#                 ${samtoolsPath}/samtools view -@ ${nthreads} ${samtoolsViewPara} - -o ${id}.bam
#             ${samtoolsPath}/samtools fixmate -@ ${nthreads} ${id}.bam ${id}_fixmate.bam
#             # rm ${id}.bam
#             ${samtoolsPath}/samtools index -@ ${nthreads} ${id}.bam

#         # bam to bed (Illumina read ID: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  <read>:<is filtered>:<control number>:<index>); 
#             #              only keep <tile>:<x-pos>:<y-pos> of the read ID
#             # columns in output: 1. chr, 2. start(read1), 3.end(read2), 4.trimmed readID (shared by R1,R2, and UMI), 5. strand (read1)
#             echo ">>>>>> Convert to bed file.......`date`"
#             ${bedtools2Path}/bedtools bamtobed -i ${id}_fixmate.bam -mate1 -bedpe | \
#                 awk '{if($1==$4&&$1!="."&&$2!=0&&$5!=0){if($2>=$5&&$9=="-"&&$10=="+"){start=$5;end=$3;strand="-"}else if($2<$5&&$9=="+"&&$10=="-"){start=$2;end=$6;strand="+"};sub("^[^:]+:[^:]+:[^:]+:[^:]+:", "", $7);print $1, start, end, $7, strand}}' OFS='\t' \
#                 > ${id}.bedpe # the second awd seen not needed!

#         # add UMI to bed (only need one CPU, but ~ 140GB/155GB RAM!)
#             # UMItab: tsv file contain two columns: seqID, UMI sequences.
#             echo ">>>>>> Add UMI in bed file.......`date`"
#             cat ${id}.bedpe | awk -v UMItab=${id}.UMItab -v ID=${id} -v OFS='\t' 'BEGIN{while(getline<UMItab){vec[$1]=$2}}{if($4 in vec){print $1, $2, $3, $5, $4, vec[$4] > ID"_UMI.bed"}else{print $0 > ID"_noUMI.bed"}}'
#     fi
# done

# ####-----------------------------------------132bp---------------------------------------------------
# ####-----------------------------Step1. alignment---------------------------------------------------
# for i in `ls ${CleanDataRoot} | grep -E "R1_trim132bp.fastq.gz"`
# do
#     id=${i%_S*_L002*}
#     if [ -f ${id}.calculating ]
#     then
#         echo ">>>>>> ${id} file exist!"   
#     else
#         echo "NA" > ${id}.calculating
#         # i="si116_C5_S1_L002_R1_trim36bp.fastq.gz"
#         echo ">>>>>> ${i}    ${id}    `date`"
#         # alignment
#             echo ">>>>>> Mapping to reference genome.......`date`"
#             ${bwaPath}/bwa mem -t ${nthreads} ${referenceFile} \
#                 ${CleanDataRoot}/${i} \
#                 ${CleanDataRoot}/${i%R1_trim132bp.fastq.gz}R3_trim132bp.fastq.gz | \
#                 # ${samtoolsPath}/samtools sort -n -@ ${nthreads} | \
#                 ${samtoolsPath}/samtools view -@ ${nthreads} -bS - -o ${id}_raw.bam
#             # ${samtoolsPath}/samtools fixmate -@ ${nthreads} ${id}.bam ${id}_fixmate.bam
#             # rm ${id}.bam
#             # ${samtoolsPath}/samtools index -@ ${nthreads} ${id}.bam

#         # # bam to bed (Illumina read ID: @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>  <read>:<is filtered>:<control number>:<index>); 
#         #     #              only keep <tile>:<x-pos>:<y-pos> of the read ID
#         #     # columns in output: 1. chr, 2. start(read1), 3.end(read2), 4.trimmed readID (shared by R1,R2, and UMI), 5. strand (read1)
#         #     echo ">>>>>> Convert to bed file.......`date`"
#         #     ${bedtools2Path}/bedtools bamtobed -i ${id}_fixmate.bam -mate1 -bedpe | \
#         #         awk '{if($1==$4&&$1!="."&&$2!=0&&$5!=0){if($2>=$5&&$9=="-"&&$10=="+"){start=$5;end=$3;strand="-"}else if($2<$5&&$9=="+"&&$10=="-"){start=$2;end=$6;strand="+"};sub("^[^:]+:[^:]+:[^:]+:[^:]+:", "", $7);print $1, start, end, $7, strand}}' OFS='\t' \
#         #         > ${id}.bedpe # the second awd seen not needed!

#         # # add UMI to bed (only need one CPU, but ~ 140GB/155GB RAM!)
#         #     # UMItab: tsv file contain two columns: seqID, UMI sequences.
#         #     echo ">>>>>> Add UMI in bed file.......`date`"
#         #     cat ${id}.bedpe | awk -v UMItab=${id}.UMItab -v ID=${id} -v OFS='\t' 'BEGIN{while(getline<UMItab){vec[$1]=$2}}{if($4 in vec){print $1, $2, $3, $5, $4, vec[$4] > ID"_UMI.bed"}else{print $0 > ID"_noUMI.bed"}}'
#     fi
# done

