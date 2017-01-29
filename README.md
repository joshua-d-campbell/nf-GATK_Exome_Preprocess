# nf-GATK_Exome_Preprocess

##Description
Adapted from the GATK best practice guide to preprocess whole exome sequencing (WES) data

##Tool documentation
####Preprocessing
Mapping: http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently


Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar


Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels

Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

####Qualiy control metrics:

FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Picard: https://broadinstitute.github.io/picard/

####Plotting:
MultiQC: multiqc.info

##Preparing reference files before running workflow
See here for downloading reference and vcf files:
https://software.broadinstitute.org/gatk/guide/article?id=1213

This workflow assumes that the FASTA reference file has been preprocessed and indices have been made according to:
http://gatkforums.broadinstitute.org/wdl/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

##Input parameters

####Target and bait information for the hybrid capture protocol:
params.targets = "SureSelect_Human_All_Exon_V4_plus_UTRs.UCSC_hg19.targets.interval_list"
params.baits = "SureSelect_XT_V4_plusUTRs_baits.UCSC_hg19.capture.interval_list"

####params.gatk_jar = "/share/pkg/gatk/3.5/install/GenomeAnalysisTK.jar"
params.picard_jar = "/share/pkg/picard/2.8.0/install/lib/picard.jar"

####Reference files:

BWA Reference file in FASTA format
params.ref = "ucsc.hg19.fasta"

####Indels and DBSNP VCFs for alignment around indels and base quality score recalibration (BQSR)
params.gold_indels1 = "1000G_phase1.indels.hg19.sites.vcf"
params.gold_indels2 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
params.dbsnp = "dbsnp_138.hg19.vcf"

####Workflow parameters
Final output directory for linked files:
params.output_dir = "outputDir"

Whether or not the input file has a header:
params.infile_header = true



