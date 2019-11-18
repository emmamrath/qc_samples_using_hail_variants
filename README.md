QC samples using Hail variants
==============================

## Overview
This pipeline uses a multi-sample joint-called VCF file to carry out quality control (QC) of those samples.  
This pipeline runs Hail to produce sample-level statistics per sample, which can be used to filter samples for QC.  

There are 2 different sample QC pipelines:
* QC1 uses [Hail 0.1](https://hail.is/docs/0.1/index.html)  and is the QC criteria used in the [MGRB paper] (https://www.biorxiv.org/content/10.1101/473348v1).
* QC2 uses [Hail 0.2](https://hail.is/docs/0.2/) and is the sample variant QC criteria used for Gnomad 2.1.1, and being used in the [ISKS](http://sarcomahelp.org/articles/sarcoma-kindred-study.html) work.

## References
* Variants are the germline genetic mutation SNVs recorded in a multi-sample [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) file.  
* [Hail](https://hail.is/docs/0.2/) is an open-source library for scalable data exploration and analysis, with a particular emphasis on genomics. Hail uses Spark to enable scalability.  
* [Gnomad 2.1.1](https://gnomad.broadinstitute.org/downloads) is a resource developed by an international coalition of investigators, with the goal of aggregating and harmonizing both exome and genome sequencing data from a wide variety of large-scale sequencing projects, and making summary data available for the wider scientific community. QC2 criteria of Gnomad 2.1.1 are detailed in [Karczewski et al. 2019](https://www.biorxiv.org/content/10.1101/531210v2).
* MGRB [(Medical Genome Reference Bank)](https://sgc.garvan.org.au/initiatives) data is a whole-genome data resource of 4000 healthy elderly individuals ([manuscript 1](https://www.biorxiv.org/content/10.1101/473348v1), [manuscript 2](https://www.nature.com/articles/s41431-018-0279-z)). QC1 criteria are detailed in [Pinese et al. 2018](https://www.biorxiv.org/content/10.1101/473348v1).  
* ISKS [(International Sarcoma Kindred Study)](http://sarcomahelp.org/articles/sarcoma-kindred-study.html) data is a whole-genome data resource (manuscript in preparation).  

### QC1 citation
The Medical Genome Reference Bank: Whole genomes and phenotype of 2,570 healthy elderly.  
Mark Pinese, Paul Lacaze, Emma M. Rath, et al.  
[https://www.biorxiv.org/content/10.1101/473348v1](https://www.biorxiv.org/content/10.1101/473348v1)

### QC2 citations
Variation across 141,456 human exomes and genomes reveals the spectrum of loss-of-function intolerance across human protein-coding genes.  
Konrad J Karczewski, Laurent C Francioli, Grace Tiao, et al.  
[https://www.biorxiv.org/content/10.1101/473348v1](https://www.biorxiv.org/content/10.1101/473348v10)  

Comparison of [ISKS (International Sarcoma Kindred Study)](http://sarcomahelp.org/articles/sarcoma-kindred-study.html) and [MGRB (Medical Genome Reference Bank)](https://sgc.garvan.org.au/initiatives) (manuscript in preparation).

### Pipeline flow
* For QC1, please follow the commands in code/commands_for_QC1.txt
* For QC2, please follow the commands in code/commands_for_QC2.txt

### QC criteria

Please note that the pipelines for QC1 and QC2 produce the necessary sample-level variant statistics for QC to be performed, but do not perform QC.  

#### QC1 sample variant criteria:  

For each cohort, choose from high-confident non-rare SNPs (eg. Illumina array SNPs were used for MGRB and ISKS studies) those SNPs that are:  
* not autosomal (not X or Y)  

These locii are used for sample QC1. Sample passes QC1 when:  
* call rate > 0.98  
* depth stdev < 10  
* heterozygous locii stdev < 1  
* rHetHomVar < 2  
* (num_singletons / (num_SNPs + num_insertions + num_deletions)) < 0.001  

#### QC2 sample variant criteria:  

For each cohort, choose from high-confident non-rare SNPs (eg. Illumina array SNPs were used for MGRB and ISKS studies) those SNPs that are:  
* biallelic (have at least one sample having a non-reference allele)  
* call_rate > 0.99  
* avg VAF > 0.001  

These locii are used for sample QC2. Sample passes QC2 when:  
* call_rate >= 0.895  
* f_stat > 0.8 or f_stat < 0.2 f_stat is the inbreeding coefficient calculated using chrom X. According to Karczewski et al. 2019, f_stat > 0.8 is classified as male, f_stat < 0.2 is classified as female. Samples with intermediate f_stat are classified as ambiguous sex.

Please note QC2 also has bam-level criteria:
* sequencing depth >= 15x
* rate of chimeric reads (sequencing read mate mapped to a different chromosome) <= 0.05
* median insert size >= a certain size such that outliers will be excluded. For Gnomad 2.1.1 this criteria was median insert size >= 250 bp.
* verifyBamID Freemix (a measure of contamination) <= 0.05


