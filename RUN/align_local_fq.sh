#!/bin/bash
###############################################################################
# Ludo Pagie, 22 March, 2012, align_local_fq_LP120405.sh
#
# DESCRIPTION:
#   wrapper around bowtie2, takes as input a single fastq file and an
#   assembly
#
# OUTPUT:
#   sorted bamfile and an index (bai) file. mapped reads are filtered
#   for alignment quality score 25 or higher
#
# VERSIONS:
#   120322: - add an exit code indicating failure/succes
#           - store output in temp files, at end move temp files to proper
#             filesnames
#           - check whether outfiles exist already, if so return succes
#           - add some mapping statistics to BWT_STATS
#   120405: - save unmapped reads in a separate file, counted and sorted
#             along abundance
#           - before mapping clip adapter sequences using cutadapt;
#             first find reads with long adapter (....
#           - remove duplicate reads from inner reads
#           - split mapped reads in 0x, 1x, 2x, 3x or more 
#   120424: - align unmapped reads to the sequence of the inserted
#             transgene (if given)
#   130117: - renamed to align_local_fq.sh for inclusion in git repos
#           - changed bowtie2 version to 2.0.5
#           - added codedir as CLI arg, should pass the directory with
#             all scripts
#   130123: - added some explanatory comments
###############################################################################

if [ $# -lt 2 ]
then
  echo 1>&2 Usage: align_local_fq.sh input_file.fq.gz assembly_version codedir [model-matrix-file]
  exit 1
fi

# set some paths for executables
BOWTIE2=bowtie2
CUTADAPT=cutadapt
FASTX_REVCOM=fastx_reverse_complement
BOWTIE2_INDEXES=~/data/DAM/indexes/

# echo some versioninfo to log:
echo 'using bowtie2 version:'
echo `${BOWTIE2} --version`
echo ''
echo 'using cutadapt version:'
echo `${CUTADAPT} --version`
echo ''
echo 'using fastx_reverse_complement version:'
echo `${FASTX_REVCOM} -h | grep FASTX`
echo ''

#set some parameters
IN_FQ=$1
MIN_Q=25
#adaptor sequences used in DamID-seq
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ADPTR_LONG_5="CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAG"
#make reverse complement of adapter sequences 
#(${FASTX_REVCOM} expects fasta input, "awk 'NR > 1'" means print everything beyond line 1)
ADPTR_SHORT_3=`echo ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
ADPTR_LONG_3=`echo ">\n${ADPTR_LONG_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`

# set base filename
case $1 in
  *.fq.gz) 
    OUT_BAM=${IN_FQ%.fq.gz}_local.bam
    CAT=zcat ;;
  *.fastq.gz) 
    OUT_BAM=${IN_FQ%.fastq.gz}_local.bam
    CAT=zcat ;;
  *.fastq) 
    OUT_BAM=${IN_FQ%.fastq}_local.bam
    CAT=cat ;;
  *.fq) 
    OUT_BAM=${IN_FQ%.fq}_local.bam
    CAT=cat ;;
  *) 
    echo "inputfile (${IN_FQ}) should either be gzipped fastq"
    echo "(*.fq.gz or fastq.gz) or fastq (*.fq or *.fastq)" 
    exit 1 ;;
esac

# set additional output files based on base filename
OUTx2_BAM=${OUT_BAM%_local.bam}_local_2x.bam
OUTx3_BAM=${OUT_BAM%_local.bam}_local_3x.bam
BWT_STATS=${OUT_BAM%_local.bam}_local.bowtie_stats
UNMAPPED_INNER=${OUT_BAM%_local.bam}_unmappedInnerReads.txt.gz
UNMAPPED_EDGE=${OUT_BAM%_local.bam}_unmappedEdgeReads.txt.gz
CLIP_STATS=${OUT_BAM%_local.bam}_local.clip_stats

# set species and assembly (for read alignment)
case $2 in
  dm3)
    ASSEMBLY=dmel_r5.41 ;;
  hg18)
    ASSEMBLY=hg18 ;;
  hg19)
    ASSEMBLY=hg19 ;;
  mm9)
    ASSEMBLY=mm9 ;;
  *)
    echo "assembly $2 not know, exiting"
    exit 1 ;;
esac

################################################
# create temporary output files
################################################
TMP_BAM_INNER=`mktemp tmp_bam_inner.XXXXXXXXXX`
TMP_BAM_EDGE=`mktemp tmp_bam_edge.XXXXXXXXXX`
TMP_STATS_INNER=`mktemp tmp_stats_inner.XXXXXXXXXX`
TMP_STATS_EDGE=`mktemp tmp_stats_edge.XXXXXXXXXX`
TMP_FQ=`mktemp tmp_fq.XXXXXXXXXX`
TMP_FQ_EDGE=`mktemp tmp_fq_edge.XXXXXXXXXX`
TMP_FQ_INNER=`mktemp tmp_fq_inner.XXXXXXXXXX`

# define function that checks exit code last command
CheckExit()
{
if [ $1 -ne 0 ]; then
  echo "$2" 1>&2
#  rm -f ${TMP_BAM_INNER}
#  rm -f ${TMP_BAM_EDGE}
#  rm -f ${TMP_STATS_INNER}
#  rm -f ${TMP_STATS_EDGE}
#  rm -f ${TMP_FQ}
#  rm -f ${TMP_FQ_EDGE}
#  rm -f ${TMP_FQ_INNER}
  exit 1
fi
}

###################################################
# awk script to split bam file in 4 files:
# unmapped, 1x-mapped, 2x-mapped, 3x-or-more mapped
###################################################
AWKSCRIPT='
BEGIN {
  flag = lshift(1, 8); # bit-8 flag
}
{
  if ( ! and ( flag, $2 ) ) { # if not first alignment of this read
    # check whether previous read is unmapped or uniquely mapped
    if ( cnt == 1 ) { if (oldflag == 4) print ID[1] > OUT0; else print ID[1] >> OUT1; }  # print to unmapped or unique mapped
    if ( cnt == 2 ) { print ID[1] >> OUT2; print ID[2] >> OUT2; } # print to 2x mapped file
    if ( cnt >= 3 ) { for (i in ID) print ID[i] >> OUT3; } # print to 3x or more mapped file
    cnt = 1; delete ID; ID[cnt] = $0; oldflag = $2; # reset intermediate variables
} else {
  cnt++; ID[cnt] = $0; # increase counter for multi mapped reads
}
}
END {
  if ( cnt == 1 ) { if (oldflag == 4) print ID[1] > OUT0; else print ID[1] >> OUT1; }
  if ( cnt == 2 ) { print ID[1] >> OUT2; print ID[2] >> OUT2; }
  if ( cnt >= 3 ) { for (i in ID) print ID[i] >> OUT3; }
}
'

# check whether outfiles exists already
if [ -f "${OUT_BAM}" -a -f "${BWT_STATS}" -a -f "${BAM_OUT}".bai ]; then
  echo "${IN_FQ} is aligned already, exiting" 1>&2
  exit 0
else
# remove if only some exist
  rm -f "${OUT_BAM}" "${BWT_STATS}" "${BAM_OUT}".bai
fi

###############################################################################
# trim adapter sequences from reads, split fastq in inner and edge reads
###############################################################################
# clip long adapter, save non-adapter reads, and adapter-match stats
ADPTR_LONG_LEN=`echo "${ADPTR_LONG_5}" | wc -c`
ADPTR_LONG_LEN=`expr ${ADPTR_LONG_LEN} - 1`
LEN_THRES=`expr ${ADPTR_LONG_LEN} / 2`
#${CAT} ${IN_FQ} | 
${CUTADAPT} -g "${ADPTR_LONG_5}" -a "${ADPTR_LONG_3}" -O "${LEN_THRES}" --match-read-wildcards --discard-trimmed ${IN_FQ} -o "${TMP_FQ}" > ${CLIP_STATS}
# clip short adapter, save non-adapter and trimmed reads in separate
# files, and save adapter-match stats
cat ${TMP_FQ} | ${CUTADAPT} -g "${ADPTR_SHORT_5}" -a "${ADPTR_SHORT_3}" --match-read-wildcards - --untrimmed-output "${TMP_FQ_INNER}" -o "${TMP_FQ_EDGE}" >> ${CLIP_STATS}
# remove fastq file with reads not trimmed with long adapters
rm ${TMP_FQ}
###############################################################################
###############################################################################

#####################################################
# run bowtie on fastq files
#####################################################
#Get number of CPU Cores
CORE=`lscpu | grep 'CPU(s):' | sed -n '1p' | rev | cut -c 1`
#Run Bowtie
BOWTIE_PAR="-k 3 -p ${CORE} -t --phred33 --local -x ${BOWTIE2_INDEXES}${ASSEMBLY}"
cat ${TMP_FQ_INNER} | (${BOWTIE2} ${BOWTIE_PAR} -U - | samtools view -bS - -o ${TMP_BAM_INNER} ) 2> ${TMP_STATS_INNER}
CheckExit $? "bowtie2 failed on inner reads"
echo '1'
# # extract all unmapped reads, collect read sequences, and read count.
# # Save in separate files
# samtools view ${TMP_BAM_INNER} | awk '$2 == 4' | cut -f 10 | sort | uniq -c |\
#   sort -gr | gzip -  > ${UNMAPPED_INNER}
#   CheckExit $? "failed to collect unmapped reads from inner reads"
#   echo '2'
#################################################################################
# split bamfile in unmapped, uniquely mapped reads, reads mapped twice, and reads
# mapped 3 times or more
# awk script above
# then prepare output samfiles with sam-header, not for OUT0.sam
samtools view -H ${TMP_BAM_INNER} > OUT1.sam
samtools view -H ${TMP_BAM_INNER} > OUT2.sam
samtools view -H ${TMP_BAM_INNER} > OUT3.sam
# and run the awk script
samtools view ${TMP_BAM_INNER} | \
awk -vOUT0=OUT0.sam -vOUT1=OUT1.sam -vOUT2=OUT2.sam -vOUT3=OUT3.sam "$AWKSCRIPT"

# unmapped reads: sort, unique, count reads only
cat OUT0.sam | cut -f 10 | sort | uniq -c |\
  sort -gr | gzip -  > ${UNMAPPED_INNER}
rm -f OUT0.sam

# filter for quality score MIN_Q
samtools view -S -b -h -q ${MIN_Q} OUT1.sam -o tmp.bam
CheckExit $? "samtools filter Q failed on inner reads"
echo '3'
# sort bam file
mv -f tmp.bam ${TMP_BAM_INNER}
rm -f OUT1.sam
samtools sort ${TMP_BAM_INNER} srt
CheckExit $? "samtools sort failed on inner reads"
echo '4'
# index bam file
mv -f srt.bam ${TMP_BAM_INNER}
# samtools index ${TMP_BAM_INNER}
# CheckExit $? "samtools index failed on inner reads"
echo '5'

## TODO
# still need to do something with OUT2.sam and OUT3.sam 

cat ${TMP_FQ_EDGE} | \
  (${BOWTIE2} ${BOWTIE_PAR} -U - | \
  samtools view -bS - -o ${TMP_BAM_EDGE} ) 2> ${TMP_STATS_EDGE}
CheckExit $? "bowtie2 failed on edge reads"
echo '11'
# extract all unmapped reads, collect read sequences, and read count.
# Save in separate files
# samtools view ${TMP_BAM_EDGE} | awk '$2 == 4' | cut -f 10 | sort | uniq -c |\
#  sort -gr | gzip - > ${UNMAPPED_EDGE}
# CheckExit $? "failed to collect unmapped reads from edge reads"
# echo '22'

# prepare OUT1.sam. OUT2.same and OUT3.sam will be appended to
samtools view -H ${TMP_BAM_EDGE} > OUT1.sam
# and run the awk script
samtools view ${TMP_BAM_EDGE} | \
awk -vOUT0=OUT0.sam -vOUT1=OUT1.sam -vOUT2=OUT2.sam -vOUT3=OUT3.sam "$AWKSCRIPT"

# unmapped reads: sort, unique, count reads only
cat OUT0.sam | cut -f 10 | sort | uniq -c |\
  sort -gr | gzip -  > ${UNMAPPED_EDGE}
rm -f OUT0.sam

# filter for quality score MIN_Q
samtools view -b -h -q ${MIN_Q} ${TMP_BAM_EDGE} -o tmp.bam
CheckExit $? "samtools filter Q failed on edge reads"
echo '33'
# sort bam file
mv -f tmp.bam ${TMP_BAM_EDGE}
rm -f OUT1.sam
samtools sort ${TMP_BAM_EDGE} srt
CheckExit $? "samtools sort failed on edge reads"
echo '44'
# index bam file
mv -f srt.bam ${TMP_BAM_EDGE}
# samtools index ${TMP_BAM_EDGE}
# CheckExit $? "samtools index failed on edge reads"
echo '55'

# compress OUT2.sam and OUT3.sam to bam
samtools view -S -b OUT2.sam -o ${OUTx2_BAM}
rm -f OUT2.sam
samtools view -S -b OUT3.sam -o ${OUTx3_BAM}
rm -f OUT3.sam

# merge stat files
cat ${TMP_STATS_INNER} ${TMP_STATS_EDGE} > ${BWT_STATS}
OUT=$?
if [ ${OUT} -ne 0 ]; then
  rm -f ${OUT_BAM}
  rm -f ${OUT_BAM}.bai
  CheckExit ${OUT} "merging stat files failed"
fi
# delete tmp stats files
rm -f ${TMP_STATS_INNER} ${TMP_STATS_EDGE}
echo '6'
#####################################################
#####################################################

#####################################################
# map unmapped reads onto transgene sequence NOT!!!!!
#####################################################
# parse model matrix file; extract filename with transgene

# construct reference sequence file for bowtie
# run bowtie

#####################################################
#####################################################

#####################################################
# collect additional statistics on reads and mapping:
#####################################################
# count reads after ${Q_MIN} filtering
######################################
CNT_INNER=`samtools view -c "${TMP_BAM_INNER}"`
CheckExit $? "samtools count inner after filter Q failed"
echo '7'
CNT_EDGE=`samtools view -c "${TMP_BAM_EDGE}"`
CheckExit $? "samtools count edge after filter Q failed"
echo '8'
echo "inner nreads after Q_filter: ${CNT_INNER}" >> ${BWT_STATS}
echo "edge nreads after Q_filter: ${CNT_EDGE}" >> ${BWT_STATS}
echo '9'
CNT_ALL=`expr ${CNT_INNER} + ${CNT_EDGE}`
CheckExit $? "summing inner and edge reads failed"
echo '10'
echo "all nreads after Q_filter: ${CNT_ALL}" >> ${BWT_STATS}
# removing duplicates from inner reads
samtools rmdup -s ${TMP_BAM_INNER} tmp.bam
CheckExit $? "samtools rmdup inner reads failed"
echo '12'
######################################
# # reads lost after rmdup
######################################
CNT_INNER_RMDUP=`samtools view -c tmp.bam`
# CNT=`samtools view -h "${TMP_BAM_INNER}" | \
#      samtools view -S -h -b - -o - | \
#      samtools rmdup -s - - | \
#      samtools view -c -`
CheckExit $? "samtools unique count inner reads failed"
echo "nreads unique inner reads: ${CNT_INNER_RMDUP}" >> ${BWT_STATS}
# replace inner reads with unique reads only
mv -f tmp.bam ${TMP_BAM_INNER}
CheckExit $? "replacing inner reads with unique reads failed"
echo '13'
# removing duplicates from edge reads
samtools rmdup -s ${TMP_BAM_EDGE} tmp.bam
CheckExit $? "samtools rmdup edge reads failed"
echo '14'
CNT_EDGE_RMDUP=`samtools view -c tmp.bam`
# CNT=`samtools view -h "${TMP_BAM_INNER}" | \
#      samtools view -S -h -b - -o - | \
#      samtools rmdup -s - - | \
#      samtools view -c -`
CheckExit $? "samtools unique count edge reads failed"
echo "nreads unique edge inner reads: ${CNT_EDGE_RMDUP}" >> ${BWT_STATS}
# don't replace edge reads with unique reads only
rm -f tmp.bam
echo '15'
######################################
# count mito reads 
######################################
samtools index ${TMP_BAM_EDGE}
CNT_EDGE_MITO=`samtools view -c ${TMP_BAM_EDGE} dmel_mitochondrion_genome`
rm -f ${TMP_BAM_EDGE}.bai
CheckExit $? "count mito reads from edge reads failed"
echo "# mito reads from edge reads: ${CNT_EDGE_MITO}" >> ${BWT_STATS}
echo '16'
samtools index ${TMP_BAM_INNER}
CNT_INNER_MITO=`samtools view -c ${TMP_BAM_INNER} dmel_mitochondrion_genome`
rm -f ${TMP_BAM_INNER}.bai
CheckExit $? "count mito reads from inner reads failed"
echo "# mito reads from inner reads: ${CNT_INNER_MITO}" >> ${BWT_STATS}
echo '17'
samtools sort ${OUTx2_BAM} tmp
mv tmp.bam ${OUTx2_BAM}
samtools index ${OUTx2_BAM}
CNT_2X_MITO=`samtools view -c ${OUTx2_BAM} dmel_mitochondrion_genome`
rm -f ${OUTx2_BAM}.bai
CheckExit $? "count mito reads from 2xmapped reads failed from ${OUTx2_BAM}"
echo "# mito reads from 2xmapped reads: ${CNT_2X_MITO=}" >> ${BWT_STATS}
echo '18'
samtools sort ${OUTx3_BAM} tmp
mv tmp.bam ${OUTx3_BAM}
samtools index ${OUTx3_BAM}
CNT_3X_MITO=`samtools view -c ${OUTx3_BAM} dmel_mitochondrion_genome`
rm -f ${OUTx3_BAM}.bai
CheckExit $? "count mito reads from 3xmapped reads failed from ${OUTx3_BAM}"
echo "# mito reads from 3xmapped reads: ${CNT_3X_MITO=}" >> ${BWT_STATS}
echo '19'

# merge temp files to final filenames
# mv ${TMP_BAM} ${OUT_BAM}
samtools merge ${OUT_BAM} ${TMP_BAM_INNER} ${TMP_BAM_EDGE}
CheckExit $? "merging bam files failed"
# delete tmp bamm files
rm -f ${TMP_BAM_INNER} ${TMP_BAM_EDGE}
# index merged bam file
samtools index ${OUT_BAM}
OUT=$?
if [ ${OUT} -ne 0 ]; then
  rm -f ${OUT_BAM}
  CheckExit ${OUT} "samtools index failed"
fi

rm -f ${TMP_BAM_INNER}
rm -f ${TMP_BAM_EDGE}
rm -f ${TMP_FQ_EDGE}
rm -f ${TMP_FQ_INNER}

# exit succesful
exit 0
