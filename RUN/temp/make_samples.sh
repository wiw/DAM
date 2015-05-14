#!/bin/bash
FWD=/home/anton/backup/tinput/bowtie/fwd/*fastq
REV=/home/anton/backup/tinput/bowtie/rev/*fastq
pairedFolder=/home/anton/backup/tpoutput/bowtie/exp/alignedReads
outputFolder=/home/anton/backup/tinput/bowtie/samples
for i in $FWD; do
i_base=$(basename $i)
paired=$(echo $i_base | sed -r "s/(.+L6_)(.*)/\1paired/")
head -n 4000000 $pairedFolder/$paired/edge_F.fastq | grep "@R" > fwd
head -n 4000000 $pairedFolder/$paired/inner_F.fastq | grep "@R" >> fwd
perl -e 'open($fh, shift(@ARGV)); $h{$_}++ while(<$fh>); while(<>){exists $h{$_}?print $_.<>.<>.<> : <>.<>.<> }' fwd $i > $outputFolder/$i_base
gzip -c $outputFolder/$i_base > $outputFolder/${i_base}.gz
done
for i in $REV; do
i_base=$(basename $i)
paired=$(echo $i_base | sed -r "s/(.+L6_)(.*)/\1paired/")
head -n 4000000 $pairedFolder/$paired/edge_R.fastq | grep "@R" > rev
head -n 4000000 $pairedFolder/$paired/inner_R.fastq | grep "@R" >> rev
perl -e 'open($fh, shift(@ARGV)); $h{$_}++ while(<$fh>); while(<>){exists $h{$_}?print $_.<>.<>.<> : <>.<>.<> }' rev $i > $outputFolder/$i_base
gzip -c $outputFolder/$i_base > $outputFolder/${i_base}.gz
done
