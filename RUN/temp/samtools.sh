flags="83 99 147 163"
GATC_GFF="/home/anton/data/DAM/COR/DmelGATCfragments-r5_LP120507.gff"
for i in flags; do
	if (( i == 83 )); then
		opt="-h -f $i"
	else
		opt="-f $i"
	fi
	samtools view $opt tmp_bam_inner.4w0DZix2mr >> inner.sam
done
samtools view -bS inner.sam > inner.bam
samtools sort -n -@ 8 inner.bam inner-sort
rm -f inner.bam inner.sam
samtools view -h inner-sort.bam | htseq-count -i ID -m intersection-strict -s no -o out.sam - "${GATC_GFF}" | gzip -c > inner.txt.gz