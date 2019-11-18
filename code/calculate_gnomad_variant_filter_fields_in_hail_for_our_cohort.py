#!/usr/bin/python
# python calculate_gnomad_variant_filter_fields_in_hail_for_our_cohort.py -infile aa -outfile aa -locii aa

# python calculate_gnomad_variant_filter_fields_in_hail_for_our_cohort.py -infiles /nvme/emmrat/new_qc/mgrb/mgrb_qc_vcf_files.txt -outfile /nvme/emmrat/new_qc/mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.tsv.bgz -locii /nvme/emmrat/new_qc/code/qc_loci.InfiniumQCArray-24v1-0_A3.bed -locii_mt_file /nvme/emmrat/new_qc/mgrb/MGRB.phase2qc.hc.norm.InfiniumQCArray_24v1_0_A3.mt -vcf_or_mt vcf

# python calculate_gnomad_variant_filter_fields_in_hail_for_our_cohort.py -infiles /nvme/emmrat/new_qc/isks/isks_qc_vcf_files.txt -outfile /nvme/emmrat/new_qc/isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.tsv.bgz -locii /nvme/emmrat/new_qc/code/qc_loci.InfiniumQCArray-24v1-0_A3.bed -locii_mt_file /nvme/emmrat/new_qc/isks/isks.qc.norm.renamed_sampleids.InfiniumQCArray_24v1_0_A3.mt -vcf_or_mt vcf

########## Hail program to get biallelic status and call rate per variant for Illumina SNPs in ISKS.
# * Import Illumina QC SNPs. Subset all ISKS variants to just the Illumina subset. Carry out the following in Hail.
# 	- Make sure they are biallelic (that at least 2 different alleles are observed).
# 	- Calculate call rate per variant for ISKS. Confirm that every Illumina variant has high call_rate > 0.99.
# 	- Calculate allele frequency. Confirm that these are common SNVs with allele frequency > 0.1%.

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
def count_different_alleles_per_variant( hl, cohort_mt ):

	# Make sure they are biallelic (that at least 2 different alleles are observed).

	annotation_expr = {}
	annotation_expr['mean_num_alt_alleles'] = hl.agg.mean(cohort_mt.GT.n_alt_alleles())
	cohort_mt = cohort_mt.annotate_rows(**annotation_expr)

	return hl, cohort_mt

######################################################
def calculate_call_rate_per_variant( hl, cohort_mt ):

	# Calculate call rate per variant for ISKS. Confirm that every Illumina variant has high call_rate > 0.99.

	cohort_mt = cohort_mt.annotate_rows(call_rate = hl.agg.fraction(hl.is_defined(cohort_mt.GT)))

	return hl, cohort_mt

######################################################
def calculate_allele_frequency_per_variant( hl, cohort_mt ):

	# Calculate allele frequency. Confirm that these are common SNVs with allele frequency > 0.1%.

	annotation_expr = {}
	annotation_expr['variant_vaf'] = hl.agg.mean(cohort_mt.GT.n_alt_alleles()) / 2
	cohort_mt = cohort_mt.annotate_rows(**annotation_expr)

	return hl, cohort_mt

######################################################
def output_metrics( cohort_mt, outfile ):

	# write out the output file of the metric calculated in this program, one line per sample

	rows_table = cohort_mt.rows()
	fields_to_drop = ['rsid', 'qual', 'filters', 'info']
	rows_table = rows_table.drop(*fields_to_drop)
	rows_table.export( outfile )

	return

######################################################
def main():

	# what input arguments have been supplied to the program

	parser = argparse.ArgumentParser(description='Read in a vcf. Load into hail. Calculate statistics in hail that will be used for variant QC. Output statistics to tab-delimited file.')
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
	#outfile = '/nvme/emmrat/new_qc/results/MGRB_ISKS_RISC16.SNPs_for_QC.InfiniumQCArray_24v1_0_A3.tsv.bgz'
	#vcf_or_mt = 'mt'

	import hail as hl
	hl.init()

	# read in the variants for samples

	if ((vcf_or_mt == "VCF") or (vcf_or_mt == "vcf")):
		hl, locii_mt = read_vcf_data_from_vcf( hl, infiles, locii, locii_mt_file )
	else:
		hl, locii_mt = read_vcf_data_from_mt( hl, locii, locii_mt_file )

	# count number of different alleles per variant

	hl, locii_mt = count_different_alleles_per_variant( hl, locii_mt )

	# calculate variant call rate per variant

	hl, locii_mt = calculate_call_rate_per_variant( hl, locii_mt )

	# calculate variant allele frequency per variant

	hl, locii_mt = calculate_allele_frequency_per_variant( hl, locii_mt )

	# output these calculated metrics for each sample

	output_metrics( locii_mt, outfile )

	locii_mt.count()

if __name__=='__main__':
    main()

