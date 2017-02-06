# nf-GATK_Exome_Preprocess

## 1. Description
Adapted from the GATK best practice guide to preprocess whole exome sequencing (WES) data. To run, use the following command:

```
nextflow nf-GATK_Exome_Preprocess.nf -c nextflow.config -with-trace
```

If the pipeline fails at any point and you fix the issue, the pipeline can be restarted with job avoidance using this command:
```
nextflow nf-GATK_Exome_Preprocess.nf -c nextflow.config -with-trace -resume
```

For more information about Nextflow commands, please see the following link:

https://www.nextflow.io/


## 2. Documentation of tools and approaches used for preprocessing
#### Preprocessing:
- Mapping: http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
- Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
- Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
- Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr

#### Qualiy control metrics:
- FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- Picard: https://broadinstitute.github.io/picard/

#### Plotting:
- MultiQC: [multiqc.info](multiqc.info) (currently using v0.9)

## 3. Preparing reference files before running workflow
See here for downloading reference and vcf files:

https://software.broadinstitute.org/gatk/guide/article?id=1213

This workflow assumes that the FASTA reference file has been preprocessed and indices have been made according to:

http://gatkforums.broadinstitute.org/wdl/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk

## 4. Input parameters
```
params.infile = "fastq_input.txt"
params.targets = "SureSelect_Human_All_Exon_V4_plus_UTRs.UCSC_hg19.targets.interval_list"
params.baits = "SureSelect_XT_V4_plusUTRs_baits.UCSC_hg19.capture.interval_list"
params.gatk_jar = "/share/pkg/gatk/3.5/install/GenomeAnalysisTK.jar"
params.picard_jar = "/share/pkg/picard/2.8.0/install/lib/picard.jar"
params.ref = "ucsc.hg19.fasta"
params.gold_indels1 = "1000G_phase1.indels.hg19.sites.vcf"
params.gold_indels2 = "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
params.dbsnp = "dbsnp_138.hg19.vcf"
params.output_dir = "outputDir"
params.infile_header = true
```

#### Input file description
The input file should be tab delimited and contain 11 columns: 
  1. INDIVIDUAl_ID - The ID of the individual from which the sample was derived.
  2. SAMPLE_ID - The ID of the sample. Note that more than one sample can come from the same individual (e.g. tumor/normal pair)
  3. LIBRARY_ID - The ID of the DNA library. Multiple sequencing libraries can be prepared from the same sample.
  4. RG_ID - Read group ID
  5. PLATFORM_UNIT - Generally is the read group ID plus the library ID
  6. PLATFORM - Sequencer (e.g. illumina)
  7. PLATFORM_MODEL - Sequencer (e.g. HiSeq2500)
  8. RUN_DATE - Date of sequencing run
  9. CENTER - Location of sequencing run
  10. R1 - Full path to Fastq file 1
  11. R2 - Full path to Fastq file 2

For more information on how to properly form read group IDs, please see the following:

https://software.broadinstitute.org/gatk/guide/article?id=6472


#### Target and bait information for the hybrid capture protocol:
targets: The genomic location of gene/exon boundaries that are supposed to be captured (in interval list format)

baits: The genomic location of the actual baits/probes used for capture (in interval list format)

#### JAR files
gatk_jar: JAR of GATK (tested with v3.5)

picard_jar: JAR of Picard Tools (tested with v2.8.0)

#### Reference files:
ref: BWA Reference file in FASTA format

gold_indels1 = High quality indel vcf file for realignment around indels

gold_indels2 = High quality indel vcf file for realignment around indels

dbsnp = dbSNP vcf used in base quality score recalibration (BQSR)

#### Other workflow parameters
output_dir: Final output directory for linked files

infile_header: Whether or not the input file has a header

## 4. Config file
The config file "nextflow.config" is included which contains all of the input paramters. To run on a cluster, you may need to change the "executor" and ".clusterOptions" for each subtask to work on your own system. If you want to change the number of cpus or memory requirements for a subtask, you will need to change the code in the main script as these requirements are currenly hard coded in the actual Linux command. To adapt NextFlow workflows to your own cluster setup, see the following link: 

https://www.nextflow.io/docs/latest/executor.html

## 5. Program versions and dependencies
This pipeline has been successfully run with the following versions
  - Picard tools v2.8.0 (requires java v1.8)
  - GATK v3.5
  - BWA v0.7.12
  - FastQC v0.11.3
  - multiqc v0.9 (using python 2.7.12)

**Important note:** These programs are currently loaded using the "module load" command. However, this will vary from system to system depending on your local setup. Therefore you may need to delete these commands and make sure these programs are accessible in your path.

## 6. File cleanup
This workflow does not currently delete the intermediate bams produced during the various steps. There after workflow completion, the follow commands will delete all unnecessary ".bam" and ".bai" files.

```
rm -fv work/*/*/*.unaligned.bam
rm -fv work/*/*/*.adapters_marked.bam
rm -fv work/*/*/*.aligned.bam
rm -fv work/*/*/*.aligned.bai
rm -fv work/*/*/*.dedup.bam
rm -fv work/*/*/*.dedup.bai
rm -fv work/*/*/*.realign.bam
rm -fv work/*/*/*.realign.bai
```





