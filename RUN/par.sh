#!/bin/bash
fq=/home/anton/backup/tinput/bowtie/*fastq.gz
log=( $fq )
count=$(echo ${#log[@]})
core=8
ADPTR_SHORT_5="GGTCGCGGCCGAG"
ILLUMINA_5="GCTCTTCCGATCT"
FASTX_REVCOM=fastx_reverse_complement
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
ILLUMINA_3=`echo -e ">\n${ILLUMINA_5}" | ${FASTX_REVCOM} | awk 'NR > 1'` #AGATCGGAAGAGC
output=/home/anton/backup/tpoutput/bowtie/paralell
for (( i=1; i<$count; i++ )); do
	(
	pos=$((i-1))
	file=$(echo ${log[$pos]})
	fname=$(basename $file)
	sname=${fname%.fastq.gz}
	if [ -d "$output/$sname" ]; then echo "directory for output exists already: $output/$sname"; else mkdir "$output/$sname"; fi
	cutadapt -g "${ADPTR_SHORT_5}" -g "${ILLUMINA_5}" -a "${ADPTR_SHORT_3}" -a "${ILLUMINA_3}" ${file} -o $output/$sname/out.fastq > $output/$sname/clip.stats
	) &
	if (( i % core == 0 )); then wait; fi
	done
	wait
