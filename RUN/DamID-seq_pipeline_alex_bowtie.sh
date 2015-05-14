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
	echo "usage: DamID-seq_pipeline_alex.sh parameterfile"
}

# set code base
CODEDIR='/home/anton/data/DAM/RUN'
ALIGN_SCRIPT="${CODEDIR}/align_local_fq_alex_bowtie_mod.sh"
READS2GATC_SCRIPT="${CODEDIR}/reads2GATC_alex_bowtie_mod.sh"
DSCR=~/data/DAM/RUN/damid_description.csv


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
# FASTDIR
if [ -z "${FAST_DIR+'xxx'}" ]
then
echo "variable FAST_DIR not set, exiting"
exit 1
fi
if [ -d "${FAST_DIR}" ]
then
echo "directory for output exists already: ${FAST_DIR}."
echo "Attempting to continue where previous analysis broke down."
else
mkdir "${FAST_DIR}"
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
	# TMP_DIR="${OUTPUT_DIR}/alignedReads_tmp"
	# mkdir ${TMP_DIR}
for df in ${FASTQ_FILES}; do
	D=`dirname "${df}"`
	B=`basename "${df}"`
	A="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"
	ln -s ${A} ${FAST_DIR}
done  
cd ${FAST_DIR}
mkdir "${OUTPUT_DIR}/alignedReads"

FASTQ_ARR=( $FASTQ_FILES )
TOTAL=$(echo ${#FASTQ_ARR[@]})
CORE=$(lscpu | grep 'CPU(s):' | sed -n '1p' | rev | cut -c 1)
for (( i=1; i<$TOTAL; i++ )); do
	(
	pos=$((i-1))
	fq=$(echo ${FASTQ_ARR[$pos]})
	fq_local=`basename ${fq}`
	fq_base=${fq_local%.fastq.gz}
	fq_human=`grep -w $fq_local $DSCR | sed -r "s/[a-zA-Z0-9_-.]*$//;s/\t//"` # human-readable name in variable, get by parse $DSCR
	fq_tag=`echo $fq_human | sed -r "s/(.*)\.(.*)\.(.*)\.(.*)/\3/;s/.*_(F|R)/\1/"` # Tag to match forward or reverse reads
	if [ "$fq_tag" = "F" ]; then
		fq_prefix=`echo $fq_human | sed -r "s/(.+)${fq_tag}.+/\1/"` # Get the name of reads from one replicate
		fq_rep=`echo $fq_human | sed -r "s/(.+)${fq_tag}(.+)/\2/"` # Get the replicate number
		fq_local=`grep -P "${fq_prefix}(F|R)${fq_rep}" $DSCR | sed -r "s/.+[[:space:]]+(.+)/\1/" | sed ':a;N;$!ba;s/\n/ /g'` # Generate list of files who involved in alinging
		fq_base=`echo $fq_base | sed -r "s/(.*?).R.{9}/\1_paired/"`
	elif [ "$fq_tag" = "R" ]; then
		echo "Skip reverse reads from $fq_human."
		continue
	fi

	bash -x ${ALIGN_SCRIPT} ${fq_local} ${ASSEMBLY} ${FAST_DIR}
	mv "${FAST_DIR}/$fq_base" "${OUTPUT_DIR}/alignedReads/"
	
	if [ $? -ne 0 ]; then
		echo "alignscript failed on fastq file; ${fq_local}, aborting" 1>&2
		exit 1
	fi
	#rm -f ${fq_local}
) &
if (( i % CORE == 0 )); then wait; fi 
done
wait

mv "${FAST_DIR}/cutadapt_statistics.csv" "${OUTPUT_DIR}/cutadapt_statistics.csv"
rm -f ${FAST_DIR}/*
cd $ORG_DIR
# move temp output dir to $OUT_DIR
# mv "${TMP_DIR}" "${OUTPUT_DIR}/alignedReads"
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
	T_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp"
	mkdir ${T_DIR}
	mkdir ${OUTPUT_DIR}/gatcReadCounts
	for fq in ${FASTQ_FILES}; do
		fq_base=`basename $fq`
		fq_human=`grep -w $fq_base $DSCR | sed -r "s/[a-zA-Z0-9_-.]*$//;s/\t//"` # human-readable name in variable, get by parse $DSCR
		fq_tag=`echo $fq_human | sed -r "s/(.*)\.(.*)\.(.*)\.(.*)/\3/;s/.*_(F|R)/\1/"` # Tag to match forward or reverse reads
		if [ "$fq_tag" = "F" ]; then
			fq_base=`echo $fq_base | sed -r "s/(.*?).R.{9}(\.fastq\.gz)/\1_paired/"`
			INFILES=${OUTPUT_DIR}/alignedReads/$fq_base/*local.bam
			# tempdir for GATC_mapping output
			TMP_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp/$fq_base"
			if [ ! -d "${TMP_DIR}" ]; then
				mkdir ${TMP_DIR}
			fi
			cd ${TMP_DIR}
			bash -x ${READS2GATC_SCRIPT} $INFILES
			if [ $? -ne 0 ]; then
				echo "GATC coverage script failed on ${f}, aborting"
				exit 1
			fi
		elif [ "$fq_tag" = "R" ]; then
			echo "Skip reverse reads from $fq_human"
			continue
		else
			fq_base=${fq_base%.fastq.gz}
			INFILES=${OUTPUT_DIR}/alignedReads/$fq_base/*local.bam
			# tempdir for GATC_mapping output
			TMP_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp/$fq_base"
			if [ ! -d "${TMP_DIR}" ]; then
				mkdir ${TMP_DIR}
			fi
			cd ${TMP_DIR}
			for f in ${INFILES}; do
				bash -x ${READS2GATC_SCRIPT} $f
				if [ $? -ne 0 ]; then
					echo "GATC coverage script failed on ${f}, aborting"
					exit 1
				fi
			done
		fi

		# INFILES=${OUTPUT_DIR}/alignedReads/$fq_base/*local.bam
		# # tempdir for GATC_mapping output
		# TMP_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp/$fq_base"
		# if [ ! -d "${TMP_DIR}" ]; then
		# 	mkdir ${TMP_DIR}
		# fi
		# cd ${TMP_DIR}
		# for f in ${INFILES}; do
		# 	bash -x ${READS2GATC_SCRIPT} $f ${CODEDIR}
		# 	if [ $? -ne 0 ]; then
		# 		echo "GATC coverage script failed on ${f}, aborting"
		# 		exit 1
		# 	fi
		# done

		cd ${ORG_DIR}
		mv -f ${TMP_DIR} ${OUTPUT_DIR}/gatcReadCounts/$fq_base
		if [ $? -ne 0 ]; then
			echo "cannot move temp output dir to ${OUTPUT_DIR}/gatcReadCounts/$fq_base, aborting" 1>&2
			exit 1
		fi
	done
fi
echo "end of Rscript for read coverage per GATC fragment"
echo ""
## end GATC coverage ########
#############################
rm -R ${T_DIR}
