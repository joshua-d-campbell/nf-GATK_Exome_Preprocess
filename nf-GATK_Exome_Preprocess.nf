#!/usr/bin/env nextflow


//############################################################################################################################
//
// Josh Campbell
// 7/20/2016
// Peforms alignment and preprocessing of paired-end exome sequencing data.
// For all samples derived from the same individual, an indel co-cleaning step will be performed on all bams jointly
//
// GATK tutorials this pipeline is based on:
// Mapping: http://gatkforums.broadinstitute.org/gatk/discussion/6483/how-to-map-and-clean-up-short-read-sequence-data-efficiently
// Marking Duplicates: http://gatkforums.broadinstitute.org/gatk/discussion/6747/how-to-mark-duplicates-with-markduplicates-or-markduplicateswithmatecigar
// Realignment around indels: http://gatkforums.broadinstitute.org/gatk/discussion/2800/howto-perform-local-realignment-around-indels
// Base Quality Score Recalibration: http://gatkforums.broadinstitute.org/gatk/discussion/2801/howto-recalibrate-base-quality-scores-run-bqsr
// 
//############################################################################################################################


// Set up global variables for requried parameters:
inputFile = file(params.infile)
inputFileHeader = params.infile_header

// Set up global variables for parameters with preset defaults:
REF = file(params.ref)
GATK = file(params.gatk_jar)
PICARD = file(params.picard_jar)
GOLD1 = file(params.gold_indels1)
GOLD2 = file(params.gold_indels2)
TARGETS = file(params.targets)
BAITS = file(params.baits)
OUTDIR = file(params.output_dir)
DBSNP = file(params.dbsnp)

logParams(params, "nextflow_parameters.txt")

VERSION = "1.2"

// Header log info
log.info "========================================="
log.info "GATK Best Practices for Exome-Seq Preprocessing v${VERSION}"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="



//#############################################################################################################
//#############################################################################################################
//
// Main
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Send FASTQ files to two processes from input file: FastQC and FastqToSam
//
// ------------------------------------------------------------------------------------------------------------
Channel.from(inputFile)
       .splitCsv(sep: '\t', header: inputFileHeader)
       .into { readPairsFastQC; readPairsFastqToSam }


// ------------------------------------------------------------------------------------------------------------
//
// Preprocess reads
// 1) Convert to BAM
// 2) Mark Adapters
//
// ------------------------------------------------------------------------------------------------------------

process runFastqToSam {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastqToSam/"
    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastqToSam
    
    output:
    set indivID, sampleID, libraryID, rgID, file(outfile) into runFastqToSamOutput

    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".unaligned.bam"
    
    """
	java -Xmx5G -XX:ParallelGCThreads=1 -jar ${PICARD} FastqToSam \
		FASTQ=${fastqR1} \
		FASTQ2=${fastqR2} \
		OUTPUT=${outfile} \
		READ_GROUP_NAME=${rgID} \
		SAMPLE_NAME=${sampleID} \
		LIBRARY_NAME=${libraryID} \
		PLATFORM_UNIT=${platform_unit} \
		PLATFORM=${platform} \
		PLATFORM_MODEL=${platform_model} \
		SEQUENCING_CENTER=${center} \
		RUN_DATE=${run_date} \
		TMP_DIR=tmp
    """
}

process runMarkIlluminaAdapters {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/MarkIlluminaAdapters/"
    
    input:
    set indivID, sampleID, libraryID, rgID, ubam from runFastqToSamOutput
    
    output:
    set indivID, sampleID, libraryID, rgID, ubam, file(outfile_bam), file(outfile_metrics) into runMarkIlluminaAdaptersOutput
	
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".adapters_marked.bam"
    outfile_metrics = sampleID + "_" + libraryID + "_" + rgID + "_adapters_metrics.txt"
            
    """
	java -Xmx5G -XX:ParallelGCThreads=1 -jar ${PICARD} MarkIlluminaAdapters \
		I=${ubam} \
		O=${outfile_bam} \
		M=${outfile_metrics} \
		TMP_DIR=tmp
    """
}




// ------------------------------------------------------------------------------------------------------------
//
// Run BWA to align reads to genome
//
// ------------------------------------------------------------------------------------------------------------

process runBWA {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/"
	
    input:
    set indivID, sampleID, libraryID, rgID, ubam, ubamxt, metrics from runMarkIlluminaAdaptersOutput
    
    output:
    set indivID, sampleID, file(outfile_bam) into runBWAOutput
    
    script:
    outfile_bam = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"
	
    """
	set -o pipefail
	java -Dsamjdk.buffer_size=131072 -Dsamjdk.compression_level=1 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:ParallelGCThreads=1 -Xmx5G -jar ${PICARD} SamToFastq \
		I=${ubamxt} \
		FASTQ=/dev/stdout \
		CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
		TMP_DIR=tmp | \
	bwa mem -M -t 14 -p ${REF} /dev/stdin | \
	java -XX:ParallelGCThreads=1 -Xmx5G -jar ${PICARD} MergeBamAlignment \
		ALIGNED_BAM=/dev/stdin \
		UNMAPPED_BAM=${ubamxt} \
		OUTPUT=${outfile_bam} \
		R=${REF} CREATE_INDEX=true ADD_MATE_CIGAR=true \
		CLIP_ADAPTERS=false \
		CLIP_OVERLAPPING_READS=true \
		INCLUDE_SECONDARY_ALIGNMENTS=true \
		MAX_INSERTIONS_OR_DELETIONS=-1 \
		PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
		ATTRIBUTES_TO_RETAIN=XS \
		TMP_DIR=tmp
		
	rm -rf tmp	
	"""	
}




// ------------------------------------------------------------------------------------------------------------
//
// Combined libraries from the same Individual/Sample to send to MarkDuplicates
//
// ------------------------------------------------------------------------------------------------------------

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])




// ------------------------------------------------------------------------------------------------------------
//
// Run Picard MarkDuplicates
// This is used to merge different libraries from the same sample or the same library run on different lanes
// Requires a lot of memory
// Need to set "ParallelGCThreads" otherwise it will "grab" extra available threads without asking (and potentially be terminated by SGE)
//
// ------------------------------------------------------------------------------------------------------------

process runMarkDuplicates {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates"
	
    input:
    set indivID, sampleID, aligned_bam_list from runBWAOutput_grouped_by_sample
    
    output:
    set indivID, sampleID, file(outfile_bam), file(outfile_bai) into runMarkDuplicatesOutput
    set indivID, file(outfile_metrics) into runMarkDuplicatesOutput_QC
    
    script:
    outfile_bam = sampleID + ".dedup.bam"
    outfile_bai = sampleID + ".dedup.bai"
    outfile_metrics = sampleID + "_duplicate_metrics.txt"	
	        
    """
	java -Xmx30g -XX:ParallelGCThreads=5 -Djava.io.tmpdir=tmp/ -jar ${PICARD} MarkDuplicates \
		INPUT=${aligned_bam_list.join(" INPUT=")} \
		OUTPUT=${outfile_bam} \
		METRICS_FILE=${outfile_metrics} \
		CREATE_INDEX=true \
		TMP_DIR=tmp
	"""  
}




// ------------------------------------------------------------------------------------------------------------
//
// Combine samples from the same Individual (e.g. tumor/normal pair) to send to runRealignerTargetCreator
//
// ------------------------------------------------------------------------------------------------------------

runMarkDuplicatesOutput_grouped_by_sample = runMarkDuplicatesOutput.groupTuple(by: [0])



// ------------------------------------------------------------------------------------------------------------
//
// Perform realignment around indels
// 1) Identify regions for realignement
// 2) Perform realignment
//
// ------------------------------------------------------------------------------------------------------------

process runRealignerTargetCreator {
    tag "${indivID}"
    publishDir "${OUTDIR}/${indivID}/Processing/RealignerTargetCreator/"
    
    input:
    set indivID, sampleID, dedup_bam_list from runMarkDuplicatesOutput_grouped_by_sample
    
    output:
    set indivID, dedup_bam_list, file(target_file) into runRealignerTargetCreatorOutput
 	
    script:
    target_file = indivID + "_target_intervals.list"
	        
    """
	java -Xmx25G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T RealignerTargetCreator \
		-R ${REF} \
		-I ${dedup_bam_list.join(" -I ")} \
		-known ${GOLD1} \
		-known ${GOLD2} \
		-o ${target_file}
	"""  
}

process runIndelRealigner {
    tag "${indivID}"
	publishDir "${OUTDIR}/${indivID}/Processing/IndelRealigner/"
	    
    input:
    set indivID, dedup_bam_list, target_file from runRealignerTargetCreatorOutput
 	    
    output:
    set indivID, file('*.realign.bam') into runIndelRealignerOutput mode flatten
    set indivID, file('*.realign.bai') into runIndelRealignerBAIOutput 
    
    script:
            
    """
	java -Xmx25G -XX:ParallelGCThreads=2 -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T IndelRealigner \
		-R ${REF} \
		-I ${dedup_bam_list.join(" -I ")} \
		-targetIntervals ${target_file} \
		-known ${GOLD1} \
		-known ${GOLD2} \
                -maxReads 500000 \
                --maxReadsInMemory 500000 \
		-nWayOut ".realign.bam"		
	"""  
}



// ------------------------------------------------------------------------------------------------------------
//
// Perform base quality score recalibration (BQSR) including
// 1) Generate a recalibration table
// 2) Generate a new table after applying recalibration
// 3) Compare differences between recalibration tables
// 4) Apply recalibration
//
// ------------------------------------------------------------------------------------------------------------

// First we need to recapture the SampleID from the filename
runIndelRealignerOutput_split = runIndelRealignerOutput.map { indivID, file -> tuple(indivID, file.baseName.replaceAll(".dedup.realign", ""), file) }

process runBaseRecalibrator {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibrator/"
	    
    input:
    set indivID, sampleID, realign_bam from runIndelRealignerOutput_split
    
    output:
    set indivID, sampleID, realign_bam, file(recal_table) into runBaseRecalibratorOutput
    
    script:
    recal_table = sampleID + "_recal_table.txt" 
       
    """
	java -XX:ParallelGCThreads=2 -Xmx30g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-o ${recal_table}
	"""
}

process runPrintReads {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

    output:
    set indivID, sampleID, file(outfile_bam), file(outfile_bai), file(outfile_bai2) into runPrintReadsOutput_for_DepthOfCoverage, runPrintReadsOutput_for_HC_Metrics, runPrintReadsOutput_for_Multiple_Metrics, runPrintReadsOutput_for_OxoG_Metrics
    set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
            
    script:
    outfile_bam = sampleID + ".clean.bam"
    outfile_bai = sampleID + ".clean.bai"
    outfile_bai2 = sampleID + ".clean.bam.bai"           
    """
	java -XX:ParallelGCThreads=2 -Xmx25g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T PrintReads \
		-R ${REF} \
		-I ${realign_bam} \
		-BQSR ${recal_table} \
		-o ${outfile_bam}
    samtools index ${outfile_bam}
    """
}    

process runBaseRecalibratorPostRecal {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibratorPostRecal/"
	    
    input:
    set indivID, sampleID, realign_bam, recal_table from runPrintReadsOutput_for_PostRecal
    
    output:
    set indivID, sampleID, recal_table, file(post_recal_table) into runBaseRecalibratorPostRecalOutput_Analyze
        
    script:
    post_recal_table = sampleID + "_post_recal_table.txt" 
       
    """
	java -XX:ParallelGCThreads=2 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T BaseRecalibrator \
		-R ${REF} \
		-I ${realign_bam} \
		-knownSites ${GOLD1} \
		-knownSites ${GOLD2} \
		-knownSites ${DBSNP} \
		-BQSR ${recal_table} \
		-o ${post_recal_table}
	"""
}	

process runAnalyzeCovariates {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/AnalyzeCovariates/"
	    
    input:
    set indivID, sampleID, recal_table, post_recal_table from runBaseRecalibratorPostRecalOutput_Analyze

	output:
	set indivID, sampleID, recal_plots into runAnalyzeCovariatesOutput
	    
    script:
    recal_plots = sampleID + "_recal_plots.pdf" 

    """
	java -XX:ParallelGCThreads=2 -Xmx5g -Djava.io.tmpdir=tmp/ -jar ${GATK} \
		-T AnalyzeCovariates \
		-R ${REF} \
		-before ${recal_table} \
		-after ${post_recal_table} \
		-plots ${recal_plots}
    """
}    





// ------------------------------------------------------------------------------------------------------------
//
// Perform a several tasks to assess QC:
// 1) Depth of coverage over targets
// 2) Generate alignment stats, insert size stats, quality score distribution
// 3) Generate hybrid capture stats
// 4) Run FASTQC to assess read quality
//
// ------------------------------------------------------------------------------------------------------------

process runDepthOfCoverage {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/DepthOfCoverage"
	    
    input:
    set indivID, sampleID, bam, bai, bai2 from runPrintReadsOutput_for_DepthOfCoverage

    output:
    file("${prefix}*") into DepthOfCoverageOutput
    
    script:
    prefix = sampleID + "."
         
    """
	java -XX:ParallelGCThreads=2 -Djava.io.tmpdir=tmp/ -Xmx5g -jar ${GATK} \
		-R ${REF} \
		-T DepthOfCoverage \
		-I ${bam} \
		--omitDepthOutputAtEachBase \
		-L ${TARGETS} \
		-ct 10 -ct 20 -ct 50 -ct 100 \
		-o ${sampleID}

	"""
}	



process runCollectMultipleMetrics {
    tag "${indivID}|${sampleID}"
 	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
 	    
    input:
    set indivID, sampleID, bam, bai, bai2 from runPrintReadsOutput_for_Multiple_Metrics

    output:
    set indivID, file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

    script:       
    prefix = sampleID + "."

    """
	java -XX:ParallelGCThreads=2 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
		PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectSequencingArtifactMetrics \
		PROGRAM=CollectBaseDistributionByCycle \
		PROGRAM=CollectQualityYieldMetrics \
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${REF} \
		DB_SNP=${DBSNP} \
		INTERVALS=${BAITS} \
		ASSUME_SORTED=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	


process runHybridCaptureMetrics {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
	    
    input:
    set indivID, sampleID, bam, bai, bai2 from runPrintReadsOutput_for_HC_Metrics

	output:
	set indivID, file(outfile) into HybridCaptureMetricsOutput mode flatten

    script:       
    outfile = sampleID + ".hybrid_selection_metrics.txt"
    target_coverage = sampleID + ".target_coverage.txt"    
    """
	java -XX:ParallelGCThreads=2 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectHsMetrics \
		INPUT=${bam} \
		OUTPUT=${outfile} \
		TARGET_INTERVALS=${TARGETS} \
		BAIT_INTERVALS=${BAITS} \
		REFERENCE_SEQUENCE=${REF} \
		PER_TARGET_COVERAGE=${target_coverage} \
		TMP_DIR=tmp
	"""
}	


process runOxoGMetrics {
    tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics"
	    
    input:
    set indivID, sampleID, bam, bai, bai2 from runPrintReadsOutput_for_OxoG_Metrics

	output:
	set indivID, file(outfile) into runOxoGMetricsOutput mode flatten

    script:       
    outfile = sampleID + ".OxoG_metrics.txt"
    
    """
	java -XX:ParallelGCThreads=2 -Xmx5g -Djava.io.tmpdir=tmp/ -jar $PICARD CollectOxoGMetrics \
		INPUT=${bam} \
		OUTPUT=${outfile} \
		DB_SNP=${DBSNP} \
		INTERVALS=${BAITS} \
		REFERENCE_SEQUENCE=${REF} \
		TMP_DIR=tmp
	"""
}	



process runFastQC {
    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/FastQC/"
	    
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center, fastqR1, fastqR2 from readPairsFastQC

    output:
    set indivID, file("*.zip") into FastQCOutput mode flatten
   	file("*.html") into FastQCOutput2
   	
    script:

    """
    fastqc -t 1 -o . ${fastqR1} ${fastqR2}
    """
}



// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------
FastQCOutput_grouped_by_indiv = FastQCOutput.groupTuple(by: [0])

process runMultiQCFastq {
    tag "${indivID}"
	publishDir "${OUTDIR}/${indivID}/QC/Fastq"
	    
    input:
    set indivID, zip_files from FastQCOutput_grouped_by_indiv
    
    output:
    set file("multiqc_fastq_file_list.txt"), file("fastq_multiqc*") into runMultiQCFastqOutput
    	
    script:
     
    """
    echo -e "${zip_files.flatten().join('\n')}" > multiqc_fastq_file_list.txt
    multiqc -n fastq_multiqc --file-list multiqc_fastq_file_list.txt
    """
}


runMarkDuplicatesOutput_QC_grouped_by_indiv = runMarkDuplicatesOutput_QC.groupTuple(by: [0])

process runMultiQCLibrary {
    tag "${indivID}"
	publishDir "${OUTDIR}/${indivID}/QC/Library"
	    
    input:
    set indivID, files from runMarkDuplicatesOutput_QC_grouped_by_indiv

    output:
    set file("multiqc_library_file_list.txt"), file("library_multiqc*") into runMultiQCLibraryOutput
    	
    script:
    """
    echo -e "${files.flatten().join('\n')}" > multiqc_library_file_list.txt
    multiqc -n library_multiqc --file-list multiqc_library_file_list.txt
    """
}


CollectMultipleMetricsOutput_grouped_by_indiv = CollectMultipleMetricsOutput.groupTuple(by: [0])
HybridCaptureMetricsOutput_grouped_by_indiv = HybridCaptureMetricsOutput.groupTuple(by: [0])
runOxoGMetricsOutput_grouped_by_indiv = runOxoGMetricsOutput.groupTuple(by: [0])

process runMultiQCSample {
    tag "${indivID}"
	publishDir "${OUTDIR}/${indivID}/QC/Sample"
	    
    input:
	set indivID, metrics_files from CollectMultipleMetricsOutput_grouped_by_indiv
    set indivID, hybrid_files from HybridCaptureMetricsOutput_grouped_by_indiv
    set indivID, oxog_files from runOxoGMetricsOutput_grouped_by_indiv
        
    output:
    set file("sample_multiqc*"), file("multiqc_sample_file_list.txt") into runMultiQCSampleOutput
    	
    script:
    """
    echo -e "${metrics_files.flatten().join('\n')}" > multiqc_sample_file_list.txt
    echo -e "${hybrid_files.flatten().join('\n')}" >> multiqc_sample_file_list.txt
    echo -e "${oxog_files.flatten().join('\n')}" >> multiqc_sample_file_list.txt
            
    multiqc -n sample_multiqc --file-list multiqc_sample_file_list.txt
    """
}





workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}




//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}
