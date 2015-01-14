#!/bin/sh
if [ $# -lt 2 ]
then
  echo 1>&2 Usage: align_local_fq.sh input_file.fq.gz assembly_version codedir [model-matrix-file]
  exit 1
fi

# set some paths for executables
BOWTIE2=bowtie2
CUTADAPT=cutadapt
FASTX_REVCOM=fastx_reverse_complement

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
ADPTR_SHORT_3=`echo -e ">\n${ADPTR_SHORT_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`
ADPTR_LONG_3=`echo -e ">\n${ADPTR_LONG_5}" | ${FASTX_REVCOM} | awk 'NR > 1'`

