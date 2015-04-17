#!/bin/bash
paired=True
basef=/home/anton/backup/tpoutput/bowtie/exp15_04_13/alignedReads_tmp/DAM-1.FCC4JPEACXX_L6_paired/paired
header="DAM1 and DAM2"
PairedStat ()
{
	TR=`grep "Total Reads" $1 | sed -r "s/.*[[:space:]]([0-9]+)/\1/"`
	MP=`grep "Matched Pairs" $1 | sed -r "s/.*[[:space:]]([0-9]+)/\1/"`
	SR=`grep "Single Reads" $1 | sed -r "s/.*[[:space:]]([0-9]+)/\1/"`
	echo "
<div class=\"row\">
<div class=\"alert\">Paired match on $2 reads</div>
<div class=\"span4\">
<h4 align=\"center\">Total reads</h4>
<p align=\"center\"><script>document.write(number_format(${TR}, 0, '.', ' '))</script></p>
</div>
<div class=\"span4\">
<h4 align=\"center\">Matched Pairs</h4>
<p align=\"center\"><script>document.write(number_format(${MP}, 0, '.', ' '))</script></p>
</div>
<div class=\"span4\">
<h4 align=\"center\">Single reads</h4>
<p align=\"center\"><script>document.write(number_format(${SR}, 0, '.', ' '))</script></p>
</div>
</div>
<p>&nbsp;</p>" >> ${basef}/fq_stat_name_report.html
}
if [ "$paired" == "True" ]; then
	echo "
	<div class=\"page-header\">
	<h1>Paired match statistic <small>for ${header}</small></h1>
	</div>" >> ${basef}/fq_stat_name_report.html

	python /usr/local/bin/paired_sequence_match.py -i " " -v ${basef}/tmp_fq_inner.F.fastq ${basef}/tmp_fq_inner.R.fastq -p ${basef}/inner_F.fastq -p ${basef}/inner_R.fastq -s ${basef}/single_inner_reads.fastq > ${basef}/paired.stat
		PairedStat ${basef}/paired.stat "inner"

	python /usr/local/bin/paired_sequence_match.py -i " " -v ${basef}/tmp_fq_edge.F.fastq ${basef}/tmp_fq_edge.R.fastq -p ${basef}/tmp_paired_edge_F.fastq -p ${basef}/tmp_paired_edge_R.fastq -s single_edge_reads.fastq > ${basef}/paired.stat


	-i " " -v ${basef}/tmp_fq_edge.F.fastq ${basef}/tmp_fq_edge.R.fastq -p ${basef}/tmp_paired_edge_F.fastq -p ${basef}/tmp_paired_edge_R.fastq -s ${basef}/single_edge_reads.fastq > ${basef}/paired.stat
		PairedStat ${basef}/paired.stat  "edge"

	python /usr/local/bin/paired_sequence_match.py -i " " -v ${basef}/single_edge_reads.fastq ${basef}/single_inner_reads.fastq -p ${basef}/paired_s.e.i_F.fastq -p ${basef}/paired_s.e.i_R.fastq > ${basef}/paired.stat
		PairedStat ${basef}/paired.stat "unmatched"

	cat ${basef}/tmp_paired_edge_F.fastq ${basef}/tmp_paired_s.e.i_F.fastq > ${basef}/edge_F.fastq
	cat ${basef}/tmp_paired_edge_R.fastq ${basef}/tmp_paired_s.e.i_R.fastq > ${basef}/edge_R.fastq

	# rm -R ${basef}/paired.stat ${basef}/tmp*.fastq ${basef}/single*.fastq

fi
