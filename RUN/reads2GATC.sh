#!/bin/bash
###############################################################################
#
# Ludo Pagie, Tuesday, May 08 2012, reads2GATC_LP120508.sh
#
# DESCRIPTION:
#   simple script to map aligned reads from a bam file onto GATC fragments. The
#   bamfile is supplied as an argument to this script.  The GATC fragments are
#   in a gff file in
#   ${PROJECT_BASE}/data/LP130118_pipeline_data/DmelGATCfragments-r5_LP120507.gff
#   The mapping is done using the HTSeq-count python script by Simon
#   Anders (http://www-huber.embl.de/users/anders/HTSeq/doc/index.html)
#
# ARGUMENTS:
#   - bamfile with mapped reads
#   - basepath of project (GATC gff file is found relative to this path)
#
# VERSIONS:
#   130117: - renamed to reads2GATC.sh for inclusion in git repos
#           - added CODEDIR CLI arg
#   130121: - corrected old file path to DmelGATCfragments-r5_LP120507.gff
#             now pointing to ${PROJECT_BASE}/data/LP130118_pipeline_data/DmelGATCfragments-r5_LP120507.gff
#   130225: - some corrections in comments
###############################################################################

if [ $# -ne 2 ]; then
  echo 'usage: reads2GATC.sh mappedReads.bam codedir'
  echo 'aborting'
  exit 1
fi

BAMFILE=$1
CODEDIR=$2

# read command line arguments, put in absolute path
# bamfile with aligned reads
D=`dirname "${BAMFILE}"`
B=`basename "${BAMFILE}"`
SAM_IN="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"
# directory path of the projects dir
PROJECT_BASE="`cd \"${CODEDIR}/..\" 2>/dev/null && pwd || echo \"${CODEDIR}/..\"`"
# gff file with GATC genomic positions
GATC_GFF="${PROJECT_BASE}/COR/DmelGATCfragments-r5_AI120515.gff"

# log some bookkeeping
echo ""
echo "Starting reads2GATC_LP120508.sh, started as:"
echo "$@"
echo ""
echo "reading from samfile: ${SAM_IN}"
echo ""
echo ""
echo 'using for projectbase:'
echo "${PROJECT_BASE}"
echo ""
echo ""
echo 'running reads2GATC.sh'
echo ''
echo 'using HTSeq-count version:'
echo `htseq-count -h | grep version`
echo ''
echo 'using R version:'
echo `R --version`
echo ''
echo ''
echo "using ${GATC_GFF} as GATC.gff file"
echo ""
echo ""


# setup additional variables;
#   - sam output file
#   - counts output file
SAM_OUT=${SAM_IN%.bam}_GATCmapped.bam
COUNTS_OUT=`basename "${SAM_IN%.bam}_GATCcounts.txt.gz"`
DF_OUT=`basename "${SAM_IN%.bam}_GATCcounts.RData"`
RD_OUT=`basename "${SAM_IN%.bam}_GATCcounts_RD.RData"`
WIG_OUT=`basename "${SAM_IN%.bam}_GATCcounts.wig"`
S_WIG_OUT=`basename "${SAM_IN%.bam}_GATCcounts_sparse.wig"`

###################################################
###  END INITIALIZATION  ##########################
###################################################


###################################################
#######  START OF READS2GATC MAPPING  #############
###################################################
# run HTSeq-count with options:
# -i ID; 'GFF attribute to be used as feature ID'. i
#        The gff file used has features like 'ID=gene:DmelGATCr5chr2L00041'
# -s no; the assay is *not* strand specific
# -o out.sam; name of output file
# -; input from stdin
# ${GATC_GFF}; the gff file containing the features to which to map
samtools view -h "${SAM_IN}" | \
  htseq-count -i ID -m intersection-strict -s no -q -o out.sam - ${GATC_GFF} | \
  gzip -c > ${COUNTS_OUT}
if [ $? -ne 0 ] ; then
  echo "HTSeq-count failed on input file ${SAM_IN}"
  exit 1
fi
###################################################
#######  END OF READS2GATC MAPPING  ###############
###################################################

##########################################
# rewrite data in various formats and save
##########################################
# convert sam.out to regular bam file
TMP_SAM=`mktemp ./tmp_bam.XXXXXXXXXX`
echo "1"
samtools view -H "${SAM_IN}" > "${TMP_SAM}"
echo "2"
cat out.sam >> "${TMP_SAM}"
echo "3"
samtools view -S -b "${TMP_SAM}" -o "${SAM_OUT}"
echo "4"
rm -f "${TMP_SAM}"
rm -f out.sam
# convert counts-file to data.frame
RSCRIPT="
# import GATC gff file
######################
GATC.fname <- '${GATC_GFF}'
colClasses <- c('character','NULL','NULL','integer','integer','NULL','NULL','NULL','character')
GATC.fragments <- read.delim(file=GATC.fname, header=FALSE, stringsAsFactors=FALSE, colClasses=colClasses)
colnames(GATC.fragments) <- c('seqname','start','end','attribute')
GATC.fragments\$attribute <- gsub('ID=gene:','',GATC.fragments\$attribute)
# import read count data
########################
counts.fname <- '${COUNTS_OUT}'
con <- gzfile(description=counts.fname)
colClasses <- c('character','integer')
GATC.counts <- read.delim(file=con, colClasses=colClasses, header=FALSE)
colnames(GATC.counts) <- c('attribute','count')
GATC.counts\$attribute <- gsub('gene:','',GATC.counts\$attribute)
# remove the additional count lines from the data files, eg;
############################################################
# 461115             no_feature      0
# 461116              ambiguous 373178
# 461117          too_low_aQual      0
# 461118            not_aligned      0
# 461119   alignment_not_unique      0
GATC.counts <- GATC.counts[GATC.counts\$attribute %in% GATC.fragments\$attribute,]

# merge counts and fragments and create data frame
##################################################
reads2GATC <- GATC.fragments
idx <- match(GATC.fragments\$attribute, GATC.counts\$attribute)
reads2GATC\$count <- GATC.counts\$count[idx]
# re-organize data frame somewhat
colnames(reads2GATC)[4] <- 'ID'
reads2GATC <- reads2GATC[, c('ID', 'seqname','start', 'end', 'count')]

# save data frame in RData format
#################################
fname <- '${DF_OUT}'
save(file=fname, x=reads2GATC)
# save as wig file
##################
wig.fname <- '${WIG_OUT}'
cat(paste('track type=wiggle_0 \n', sep = ''), file = wig.fname, append = FALSE);
WriteWig <- function(counts, out.file=wig.fname, sparse=FALSE) {
  # write a data frame 'counts' containing seqname/start/end/count as wigfile
  # the 'counts' data frame should contain a single seqname
  if(length(unique(counts\$seqname)) > 1) {
    stop('WriteWig(..) called with a data frame continaing more than 1 seqname')
  }
  seqname <- counts\$seqname[1]
  cat(paste('variableStep chrom=', seqname, '\n', sep = ''), file = out.file, append = TRUE)
  if (sparse) {
    zeros <- counts\$count == 0
    zeros.diff <- diff(zeros)
    zeros.diff <- c(1, zeros.diff)
    zeros.diff[length(zeros.diff)] <- 1
    zero.stretches <- which(zeros.diff==0 & counts\$count==0)
    if(length(zero.stretches) > 0) {
      counts <- counts[-zero.stretches,]
    }
  }
  write.table(
              # Write only 2 columns of wigdf.  
              counts[,c('start','count')], 
              file = out.file, 
              sep = ' ', 
              row.names = FALSE, 
              col.names = FALSE, 
              quote = FALSE, 
              dec = '.', 
              append = TRUE)
}
invisible(tapply(seq.int(nrow(reads2GATC)), reads2GATC\$seqname, function(idx) WriteWig(reads2GATC[idx,], sparse=FALSE)))

# save as sparse wig file (ie merge fragments with zero count)
##############################################################
wig.fname <- '${S_WIG_OUT}'
cat(paste('track type=wiggle_0 \n', sep = ''), file = wig.fname, append = FALSE);
invisible(tapply(seq.int(nrow(reads2GATC)), reads2GATC\$seqname, function(idx) WriteWig(reads2GATC[idx,], sparse=TRUE)))

# convert to RangedData and save
################################
library(IRanges)
irl <- IRangesList(tapply(seq.int(nrow(reads2GATC)), reads2GATC\$seqname,
    function(idx) IRanges(start=reads2GATC\$start[idx], end=reads2GATC\$end[idx],
    names=reads2GATC\$ID[idx])))
reads2GATC <- RangedData(irl, reads2GATC\$count)
rd.fname <- '${RD_OUT}'
save(file=rd.fname, x=reads2GATC)

# done
######
quit('no')
"
#echo "${RSCRIPT}"

NOW=`date +%y%m%d%H%M`
# echo "now = ${NOW}"
echo "${RSCRIPT}" | R --vanilla > "reads2GATC_Rscript_${NOW}.Rout"

exit 0;

