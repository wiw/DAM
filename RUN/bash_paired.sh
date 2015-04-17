paired_sequence_match.py -i " " -v inner-r1.fastq inner-r2.fastq -p inner_F.fastq -p inner_R.fastq -s single_inner_reads.fastq > paired.stat
paired_sequence_match.py -i " " -v edge-r1.fastq edge-r2.fastq -p paired_edge_F.fastq -p paired_edge_R.fastq -s single_edge_reads.fastq > paired.stat
paired_sequence_match.py -i " " -v single_edge_reads.fastq single_inner_reads.fastq -p paired_s.e.i_F.fastq -p paired_s.e.i_R.fastq -s single_s.e.i_reads.fastq > paired.stat
cat paired_edge_F.fastq paired_s.e.i_F.fastq > edge_F.fastq
cat paired_edge_R.fastq paired_s.e.i_R.fastq > edge_R.fastq
