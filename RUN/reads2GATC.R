GATC_GFF <- "/home/anton/data/DAM/COR/DmelGATCfragments-r5_LP120507.gff"
COUNTS_OUT_i <- "DAM-1.FCC4JPEACXX_L6_paired_inner_local_GATCcounts.txt.gz"
DF_OUT_i <- "DAM-1.FCC4JPEACXX_L6_paired_inner_local_GATCcounts.RData"
WIG_OUT_i <- "DAM-1.FCC4JPEACXX_L6_paired_inner_local_GATCcounts_RD.RData"
S_WIG_OUT_i <- "DAM-1.FCC4JPEACXX_L6_paired_inner_local_GATCcounts.wig"
RD_OUT_i <- "DAM-1.FCC4JPEACXX_L6_paired_inner_local_GATCcounts_sparse.wig"

COUNTS_OUT_e <- "DAM-1.FCC4JPEACXX_L6_paired_edge_local_GATCcounts.txt.gz"
DF_OUT_e <- "DAM-1.FCC4JPEACXX_L6_paired_edge_local_GATCcounts.RData"
WIG_OUT_e <- "DAM-1.FCC4JPEACXX_L6_paired_edge_local_GATCcounts_RD.RData"
S_WIG_OUT_e <- "DAM-1.FCC4JPEACXX_L6_paired_edge_local_GATCcounts.wig"
RD_OUT_e <- "DAM-1.FCC4JPEACXX_L6_paired_edge_local_GATCcounts_sparse.wig"
# import GATC gff file
######################
GATC.fname <- GATC_GFF
colClasses <- c('character','NULL','NULL','integer','integer','NULL','NULL','NULL','character')
GATC.fragments <- read.delim(file=GATC.fname, header=FALSE, stringsAsFactors=FALSE, colClasses=colClasses)
colnames(GATC.fragments) <- c('seqname','start','end','attribute')
GATC.fragments$attribute <- gsub('ID=gene:','',GATC.fragments$attribute)
# import read count data
########################
counts.fname1 <- COUNTS_OUT_e
con1 <- gzfile(description=counts.fname)
colClasses1 <- c('character','integer')
GATC.counts1 <- read.delim(file=con1, colClasses=colClasses1, header=FALSE)
colnames(GATC.counts1) <- c('attribute','count')
GATC.counts1$attribute <- gsub('gene:','',GATC.counts1$attribute)
# remove the additional count lines from the data files, eg;
############################################################
# 461115             no_feature      0
# 461116              ambiguous 373178
# 461117          too_low_aQual      0
# 461118            not_aligned      0
# 461119   alignment_not_unique      0
GATC.counts1 <- GATC.counts1[GATC.counts1$attribute %in% GATC.fragments$attribute,]

# merge counts and fragments and create data frame
##################################################
reads2GATC1 <- GATC.fragments
idx1 <- match(GATC.fragments$attribute, GATC.counts1$attribute)
reads2GATC1$count <- GATC.counts1$count[idx1]
# re-organize data frame somewhat
colnames(reads2GATC1)[4] <- 'ID'
reads2GATC1 <- reads2GATC1[, c('ID', 'seqname','start', 'end', 'count')]

# save data frame in RData format
#################################
fname <- DF_OUT
save(file=fname, x=reads2GATC)
# save as wig file
##################
wig.fname <- WIG_OUT
cat(paste('track type=wiggle_0 \n', sep = ''), file = wig.fname, append = FALSE);
WriteWig <- function(counts, out.file=wig.fname, sparse=FALSE) {
  # write a data frame 'counts' containing seqname/start/end/count as wigfile
  # the 'counts' data frame should contain a single seqname
  if(length(unique(counts$seqname)) > 1) {
    stop('WriteWig(..) called with a data frame continaing more than 1 seqname')
  }
  seqname <- counts$seqname[1]
  cat(paste('variableStep chrom=', seqname, '\n', sep = ''), file = out.file, append = TRUE)
  if (sparse) {
    zeros <- counts$count == 0
    zeros.diff <- diff(zeros)
    zeros.diff <- c(1, zeros.diff)
    zeros.diff[length(zeros.diff)] <- 1
    zero.stretches <- which(zeros.diff==0 & counts$count==0)
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
invisible(tapply(seq.int(nrow(reads2GATC)), reads2GATC$seqname, function(idx) WriteWig(reads2GATC[idx,], sparse=FALSE)))

# save as sparse wig file (ie merge fragments with zero count)
##############################################################
wig.fname <- S_WIG_OUT
cat(paste('track type=wiggle_0 \n', sep = ''), file = wig.fname, append = FALSE);
invisible(tapply(seq.int(nrow(reads2GATC)), reads2GATC$seqname, function(idx) WriteWig(reads2GATC[idx,], sparse=TRUE)))

# convert to RangedData and save
################################
library(IRanges)
irl <- IRangesList(tapply(seq.int(nrow(reads2GATC)), reads2GATC$seqname,
    function(idx) IRanges(start=reads2GATC$start[idx], end=reads2GATC$end[idx],
    names=reads2GATC$ID[idx])))
reads2GATC <- RangedData(irl, reads2GATC$count)
rd.fname <- 'RD_OUT
save(file=rd.fname, x=reads2GATC)

# done
######
quit('no')