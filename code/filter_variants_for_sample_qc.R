# Rscript filter_variants_for_sample_qc.R ../mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.tsv ../mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.filtered_variants.bed
# Rscript filter_variants_for_sample_qc.R ../isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.tsv ../isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.filtered_variants.bed

options(width=240)
options(stringsAsFactors = FALSE)
library(data.table)
require(reshape)

args = commandArgs(trailingOnly=TRUE) # for production
# args=c( '../isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.tsv', '../isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.filtered_variants.tsv' )

infile = args[1]
outfile = args[2]

indata = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="" )
indata2 = transform( indata, locus2=colsplit( locus, split="\\:", names=c('a','b') ))
names(indata2)[names(indata2) == 'locus2.a'] = 'chrom'
names(indata2)[names(indata2) == 'locus2.b'] = 'pos'
indata2$pos_minus_one = as.numeric(indata2$pos) - 1
indata = indata2

indata2 = indata[ ((indata$call_rate > 0.99) & (indata$variant_vaf > 0.001) & (indata$mean_num_alt_alleles > 0)), c( "chrom", "pos_minus_one", "pos" ) ]

write.table( indata2, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )






