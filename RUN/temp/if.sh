#!/bin/sh
OUTPUT_DIR=/home/anton/data/OUT
INFILES=${OUTPUT_DIR}/alignedReads/*bam

  # tempdir for GATC_mapping output
  TMP_DIR="${OUTPUT_DIR}/gatcReadCounts"
  cd ${TMP_DIR}
  for f in ${INFILES}; do
    # do not process bam files with 2x and 3x reads (multi-reads)
    if [[ ${f} =~ "2x.bam" ]] || [[ ${f} =~ "3x.bam" ]]; then
      echo "skipping ${f}"
      continue;
    fi
done

basef=${f##*/}
# and remove extension
basef=${basef%.*}
echo "$f"
#echo "$basef"
