#!/usr/bin/python
# python calculate_gnomad_sample_filter_fields_in_hail_for_our_cohort.py -infile aa -outfile aa -locii aa

# python calculate_gnomad_sample_filter_fields_in_hail_for_our_cohort.py -outfile /nvme/emmrat/new_qc/mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.filtered_variants.sample_qc_stats.tsv.bgz -locii /nvme/emmrat/new_qc/mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.filtered_variants.bed -locii_mt_file /nvme/emmrat/new_qc/mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.mt -vcf_or_mt mt

# python calculate_gnomad_sample_filter_fields_in_hail_for_our_cohort.py -outfile /nvme/emmrat/new_qc/isks/isks.qc.norm.renamed_sampleids.filtered_InfiniumQCArray_24v1_0_A3.sample_qc_stats.tsv.bgz -locii /nvme/emmrat/new_qc/isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.filtered_variants.bed -locii_mt_file /nvme/emmrat/new_qc/isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.mt -vcf_or_mt mt

# For each sample in ISKS, for the Illumina subset of SNPs:
#    * In Hail, calculate the inbreeding coefficient (F) for these common variants on chromosome X. 
#    * Samples with F > 0.8 will be classified as male and samples with F < 0.2 will be classified as female. Samples with intermediate F values will be classified as ambiguous sex.
# For each sample in ISKS, for the Illumina subset of SNPs:
#    * In Hail, calculate sample call_rate. Call_rate < 0.895 will be considered a low call_rate.

__author__ = 'Emma M. Rath'
__copyright__ = 'Copyright 2019, Garvan Institute of Medical Research and Kinghorn Cancer Centre'

######################################################
# export PATH="/my/path/anaconda_python3_7/bin:$PATH"
# conda activate hail
# python

import sys
import argparse
sys.path.append('/nvme/emmrat/software_various/hail_python_packages')

#from gnomad_hail import *
#from gnomad_hail.resources.sample_qc import *

######################################################
def read_vcf_data_from_vcf( hl, infiles, locii, locii_mt_file ):

	# read in the vcf files

	infiles_list = []
	infiles_handle = open(infiles, 'r')
	for infile in infiles_handle:
		infile = infile.strip()
		infiles_list.append(infile)

	all_mt = hl.import_vcf(infiles_list, reference_genome='GRCh37')

	# Subset to only the locii to be used for sample QC.
	# Intervals in UCSC BED files are 0-indexed and half open. The line "5 100 105" correpsonds to the interval [5:101-5:106) in Hail's 1-indexed notation.

	interval_table = hl.import_bed(locii, reference_genome='GRCh37')
	locii_mt = all_mt.filter_rows(hl.is_defined(interval_table[all_mt.locus]))

	# Split multi-allelic variants

	locii_mt = hl.split_multi_hts(locii_mt, keep_star=False, left_aligned=True)

	# Write out the vcf data as a Hail 0.2 matrix table. Will be quicker to read in when run this program in the future.

	locii_mt.write(locii_mt_file, overwrite=True)

	return hl, locii_mt

######################################################
def read_vcf_data_from_mt( hl, locii, cohort_mt_file ):

	# read in the mt file that is the vcf variants already filtered for locii

	locii_mt = hl.read_matrix_table(cohort_mt_file)

	return hl, locii_mt

######################################################
def calculate_call_rate_per_sample( hl, cohort_mt ):

	# calculate sample call_rate. Call_rate < 0.895 will be considered a low call_rate.

	cohort_mt = cohort_mt.annotate_cols(call_rate = hl.agg.fraction(hl.is_defined(cohort_mt.GT)))

	return hl, cohort_mt

######################################################
def calculate_inbreeding_coefficient( hl, cohort_mt ):

	male_threshold = float(0.8)
	female_threshold = float(0.2)

	# Imputes sex, and annotates mt with this data.
	# Calculate the inbreeding coefficient (F) for these common variants on chromosome X. 
	# Samples with F > 0.8 will be classified as male and samples with F < 0.2 will be classified as female. Samples with intermediate F values will be classified as ambiguous sex.

	chrX_mt = hl.filter_intervals(cohort_mt, [hl.parse_locus_interval('X')])
	sex_ht = hl.impute_sex(chrX_mt.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
	sex_colnames = ['f_stat', 'is_female']
	sex_ht = sex_ht.select(*sex_colnames)
	cohort_mt = cohort_mt.annotate_cols(**sex_ht[cohort_mt.col_key])

	# These annotations will be done later because we don't have normalized_y_coverage as input.
	# cohort_mt.describe()
	# cohort_mt = cohort_mt.annotate_cols(ambiguous_sex=((cohort_mt.f_stat >= 0.5) & (hl.is_defined(cohort_mt.normalized_y_coverage) & (cohort_mt.normalized_y_coverage <= 0.1))) |
	#					(hl.is_missing(cohort_mt.f_stat)) |
	#					((cohort_mt.f_stat >= 0.4) & (cohort_mt.f_stat <= 0.6) & (hl.is_defined(cohort_mt.normalized_y_coverage) & (cohort_mt.normalized_y_coverage > 0.1))),
	#			sex_aneuploidy=(cohort_mt.f_stat < 0.4) & hl.is_defined(cohort_mt.normalized_y_coverage) & (cohort_mt.normalized_y_coverage > 0.1))

	return hl, cohort_mt

######################################################
def output_metrics( cohort_mt, outfile ):

	# write out the output file of the metric calculated in this program, one line per sample

	cols_table = cohort_mt.cols()
	# fields_to_drop = ['call_rate, 'f_stat', 'is_female']
	# cols_table = cols_table.drop(*fields_to_drop)
	cols_table.export( outfile )

	return

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in a vcf. Load into hail. Calculate statistics in hail that will be used for sample QC. Output statistics to tab-delimited file.')
	parser.add_argument('-infiles', action="store", dest="infiles", required=False, help='Input file containing list of input vcf files, one line per vcf file')
	parser.add_argument('-outfile', action="store", dest="outfile", required=True, help='Output tab-delimited file')
	parser.add_argument('-locii', action="store", dest="locii", required=True, help='Input locii, one line per locus, each line is 0-based format chr:start-end (start is 1 bp before the included start position, end is the included end position')
	parser.add_argument('-locii_mt_file', action="store", dest="locii_mt_file", required=True, help='Name of matrix table. When vcf_or_mt is vcf, then locii_mt_file will be created from input vcf. When vcf_or_mt is mt, then locii_mt_file is read in.')
	parser.add_argument('-vcf_or_mt', action="store", dest="vcf_or_mt", required=True, help='vcf when input is a file listing vcf files. mt when input the hail 0.2 matrix table already created')
	args = parser.parse_args()

	infiles = ''
	if (args.infiles is not None):
		infiles = str(args.infiles)
	locii = str(args.locii)
	locii_mt_file = str(args.locii_mt_file)
	outfile = str(args.outfile)
	vcf_or_mt = str(args.vcf_or_mt)
	#infiles = '/nvme/emmrat/new_qc/code/list_of_sample_vcfs.txt'
	#locii = '/nvme/emmrat/new_qc/test/qc_loci.InfiniumQCArray-24v1-0_A3.bed'
	#locii_mt_file = '/nvme/emmrat/new_qc/results/MGRB_ISKS_RISC16.InfiniumQCArray_24v1_0_A3.mt'
	#outfile = '/nvme/emmrat/new_qc/results/MGRB_ISKS_RISC16.WholeGenome.WholeGenome.InfiniumQCArray_24v1_0_A3.tsv.bgz'
	#vcf_or_mt = 'vcf'

	import hail as hl
	hl.init()

	# read in the variants for samples

	if ((vcf_or_mt == "VCF") or (vcf_or_mt == "vcf")):
		hl, locii_mt = read_vcf_data_from_vcf( hl, infiles, locii, locii_mt_file )
	else:
		hl, locii_mt = read_vcf_data_from_mt( hl, locii, locii_mt_file )

	# calculate variant call rate per sample

	hl, locii_mt = calculate_call_rate_per_sample( hl, locii_mt )

	# calculate calculate the inbreeding coefficient (F) for these variants on chromosome X for each sample

	hl, locii_mt = calculate_inbreeding_coefficient( hl, locii_mt )

	# output these calculated metrics for each sample

	output_metrics( locii_mt, outfile )

	locii_mt.count()

if __name__=='__main__':
    main()

