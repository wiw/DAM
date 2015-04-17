#!/bin/bash
########################################################################
# Ludo Pagie, March 6, 2012, DamID-seq_pipeline_LP120306.sh
#
# DESCRIPTION:
#   a pipeline for preprocessing fastq files from DamID-seq experiments
#   the pipeline maps reads onto GATC fragments and writes counts
#   tables, and generates several QA plots and tables
#
# DATA:
#   input data is specified in a parameter file which is supplied as
#   argument to this runner script
#
# OUTPUT:
#   data files with read counts per GATC fragment, in bed formatted
#   files and in RData files
#   several tables and plots with QA data
#
# VERSIONS:
#   120322: make the script robust against interruptions. ie:
#   - for each substage create temp output dirs, only move temp dirs to
#     final dirs if substage is finished. Now we can simply check at
#     which stage we crashed by checking tempdirs
#   - same for output files
#   - for each step check whether outfile/outdir exists in final form or
#     in temp form
#   120323: minor changes
#   120411: use new script for R analysis; Rscripts/bam2GATC_cli_LP120411.R
#     (implements just a bugfix)
#   120426: use updated align script (align_local_fq_LP120424.sh)
#     splits mapped reads in unique and multi-reads
#     collects some additional statistics
#   120508: use HTSeq-count to map reads onto GATC fragments
#   120822: trying to output entire command + args to log file
#   130116: copy to
#     zircon:/home/ludo/projects/LP130116_DamID-seq-fly-pipeline/code to
#     clean up the whole pipeline
#     added codedir arg to ALIGN_SCRIPT call
#   130121: corrected two filepaths to align_local_fq.sh and
#     reads2GATC.sh
#     - added logging of versionID of current git commit in project home
#       directory
#     
########################################################################

usage ()
{
  echo "usage: DamID-seq_pipeline_LP120321.sh parameterfile"
}

# set code base
CODEDIR='/home/anton/data/DAM/RUN'
ALIGN_SCRIPT="${CODEDIR}/align_local_fq.sh"
READS2GATC_SCRIPT="${CODEDIR}/reads2GATC.sh"

ORG_DIR=$PWD
ALLSPECIES='human fly'
ALL_ASSEMBLY_HUMAN='hg18 hg19'
ALL_ASSEMBLY_FLY='dm3 dm4'

########################################################################
######## CHECK PARAMETER FILE ##########################################
########################################################################
if [ $# -eq 0 ]
then
  echo "no parameterfile given"
  usage
  exit 1
fi
PARFILE=$1
if [ ! -f ${PARFILE} ]
then
  echo "parfile ${PARFILE} does not exist"
  usage
  exit 1
else
  echo "using par file ${PARFILE}"
fi

# redirect stdout and stderr to a log file
NOW=`date +%d-%m-%Y_%H:%M`
LOG="${PARFILE}_${NOW}.log"
exec &> ${LOG}
# echo some stats to log file
echo "running DamID-seq pipeline, with git version:"
echo `git describe` 
echo `git log --pretty=format:'%H' -n 1` 
echo "" 
echo "date = `date`" 
echo "pwd = `pwd`"
echo "commandline ="
echo "$0 $*" 
echo ""
echo ""

echo "start of initialisation"
echo ""
# import param file
. ./${PARFILE}

# check parameters
# fastq files
if [ -z "${FASTQ_FILES+'xxx'}" ]
then
  echo "variable FASTQ_FILES not set, exiting"
  exit 1
fi
for fq in ${FASTQ_FILES}
do
  if [ ! -e $fq ]
  then
    echo "file ${fq} does not exist, exiting"
    exit 1
  fi
done
# SPECIES
if [ -z "${SPECIES+'xxx'}" ]
then
  echo "variable SPECIES not set, exiting"
  exit 1
fi
correct=0
for sp in ${ALLSPECIES}
do
  echo "specs ${sp}"
  if [ "${SPECIES}" = "${sp}" ]
  then
    correct=1
  fi
done
if [ ${correct} -eq 0 ]
then
  echo "SPECIES should be in ${ALLSPECIES}, exiting"
  exit 1
fi
# ASSEMBLY
if [ -z "${ASSEMBLY+'xxx'}" ]
then
  echo "variable ASSEMBLY not set, exiting"
  exit 1
fi
if [ "${SPECIES}" = 'human' ]
then
  ALL_ASSEMBLY=${ALL_ASSEMBLY_HUMAN}
elif [ ${SPECIES} = 'fly' ]
then
  ALL_ASSEMBLY=${ALL_ASSEMBLY_FLY}
fi
for as in ${ALL_ASSEMBLY}
do
  echo "ass $as"
  if [ "${ASSEMBLY}" = "${as}" ]
  then
    correct=1
  fi
done
if [ ${correct} -eq 0 ]
then
  echo "for species ${SPECIES} ASSEMBLY should be in ${ALL_ASSEMBLY}, exiting"
  exit 1
fi
# OUTPUTDIR
if [ -z "${OUTPUT_DIR+'xxx'}" ]
then
  echo "variable OUTPUT_DIR not set, exiting"
  exit 1
fi
if [ -d "${OUTPUT_DIR}" ]
then
  echo "directory for output exists already: ${OUTPUT_DIR}."
  echo "Attempting to continue where previous analysis broke down."
else
  mkdir "${OUTPUT_DIR}"
fi

OUTPUT_DIR="`cd \"$OUTPUT_DIR\" 2>/dev/null && pwd || echo \"$OUTPUT_DIR\"`"
echo "end of initialisation"
echo ""

########################################################################
###  END OF CHECK PARAMETER FILE #######################################
########################################################################


#############################
## run bowtie ###############
echo "run bowtie2 for aligning reads to genome"
# check whether we need to restart after failed run
if [ -d ${OUTPUT_DIR}/alignedReads ]; then
  echo "bowtie has been run succesful previously, skipping this step"
else
  # tempdir for bowtie output
  TMP_DIR="${OUTPUT_DIR}/alignedReads_tmp"
  mkdir ${TMP_DIR}
  # mkdir ${OUTPUT_DIR}/alignedReads
  for df in ${FASTQ_FILES}; do
    D=`dirname "${df}"`
    B=`basename "${df}"`
    A="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"
    ln -s ${A} ${TMP_DIR}
  done  
  cd ${TMP_DIR}
#  BASE=${HOME}/projects/DamID-seq
#  SCRIPT=${BASE}/LP111111_alignedReads/align_local_fq_LP120424.sh
  for fq in ${FASTQ_FILES}; do
    fq_local=`basename ${fq}`
    bash -x ${ALIGN_SCRIPT} ${fq_local} ${ASSEMBLY} ${CODEDIR}
    if [ $? -ne 0 ]; then
      echo "alignscript failed on fastq file; ${fq_local}, aborting" 1>&2
      exit 1
    fi
    rm -f ${fq_local}
  done
  cd $ORG_DIR
  # move temp output dir to $OUT_DIR
  mv "${TMP_DIR}" "${OUTPUT_DIR}/alignedReads"
fi
echo "end of bowtie2"
echo ""
## end bowtie ###############
#############################

#############################
## run GATC coverage ########
# check whether we need to restart at later stage
echo "run Rscript for read coverage per GATC fragment"
if [ -d ${OUTPUT_DIR}/gatcReadCounts ]; then
  echo "GATC coverage has been run succesful previously, skipping this step"
else
#  BASE=${HOME}/projects/DamID-seq
#  SH_SCRIPT="${BASE}/LP120507_HTSeq_test/reads2GATC_LP120508.sh"
#  READS2GATC_SCRIPT="${CODEDIR}/reads2GATC_LP120508.sh"
 # INFILES=`ls ${OUTPUT_DIR}/alignedReads/*bam | grep -v '.*x.bam$' | grep -v '.*x_GATCmapped.bam$'`
  
  INFILES=${OUTPUT_DIR}/alignedReads/*bam

  # tempdir for GATC_mapping output
  TMP_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp"
  if [ ! -d "${TMP_DIR}" ]; then
    mkdir ${TMP_DIR}
  fi
  cd ${TMP_DIR}
  for f in ${INFILES}; do
    # do not process bam files with 2x and 3x reads (multi-reads)
     if [[ ${f} =~ "2x.bam" ]] || [[ ${f} =~ "3x.bam" ]]; then
      echo "skipping ${f}"
      continue;
     fi
    # echo "processing ${f}"
    # remove path from filename
     basef=${f##*/}
    # and remove extension
     basef=${basef%.*}
    # echo "logging to ${basef}"
    bash -x ${READS2GATC_SCRIPT} $f ${CODEDIR}
    if [ $? -ne 0 ]; then
      echo "GATC coverage script failed on ${f}, aborting"
      exit 1
    fi
  done
  cd ${ORG_DIR}
  mv -f ${TMP_DIR} ${OUTPUT_DIR}/gatcReadCounts
  if [ $? -ne 0 ]; then
    echo "cannot move temp output dir to ${OUTPUT_DIR}/gatcReadCounts, aborting" 1>&2
    exit 1
  fi
fi
echo "end of Rscript for read coverage per GATC fragment"
echo ""
## end GATC coverage ########
#############################


#############################
## process results ##########

# in R script:
# - mapping statistics (parse output of bowtie)
# - proportion duplicate reads (should be done in bowtie script; are there
#   tools other than samtools capable of collecting stats?)
# - proportion of reads spanning GATC sites (count reads containing
#   GATC (easy, but edge reads always contain GATC), or base it on
#   mapping info (maybe use first filter and then consider mapping info)
# - cumulative coverage plot of GATC fragment coverage
# - histogram; reads per fragment
# - create RangedData objects with data

# is it possible to collect these statistics over runs in some database
# so we can also show statistics of current run relative to previous
# statistics?

## process results ##########
#############################

#/usr/local/src/R-2.14.1-build_LP120202/bin/Rscript -e "library(bioinfR); parfile <- '$PARFILE'; runSweave('$1', driver='knit')"

# instead of runSweave use knit() directly plus some form of texi2pdf

