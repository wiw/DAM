CORE=`lscpu | grep 'CPU(s):' | sed -n '1p' | rev | cut -c 1`
split -a 1 -d -l $(( $( wc -l < tmp_fq_inner.F.fastq ) / $CORE )) ${fq_base} $basef/input.fastq
SPLITTER=( $(ls *.fastq?) )
for i in ${SPLITTER[@]}; do
	{
cutadapt -g "${ADPTR_SHORT_5}" -g "${ILLUMINA_5}" -a "${ADPTR_SHORT_3}" -a "${ILLUMINA_3}" -e 0.01 -n 3 --overlap 12 --match-read-wildcards --untrimmed-output $basef/untrim_out.fastq ${i} -o $basef/out.fastq > $stats/clip_${fq_base%.fastq.gz}.stats
}; &
wait
done
