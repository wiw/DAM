#!/bin/sh
INFILES=`ls /home/anton/data/OUT/alignedReads/*bam`
REGEXP_2X='2x.bam'
REGEXP_3X='3x.bam'
  for f in ${INFILES}; do
    if [[ ${f} =~ $REGEXP_2X ]] || [[ ${f} =~ $REGEXP_3X ]]; then
      echo "skipping ${f}"
      continue;
    fi
	echo "$BASH_REMATCH"
    echo $f
    done
