#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params {
    // Input/Output
    reads = "raw_data/fastq/*_L00{1,2,3,4}_R{1,2,3}_001.fastq.gz"
    outdir = "raw_data/results"

    //Reference files
    genome = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Sequence/bwa/UU_Cfam_GSD_1.0_ROSY.fa"
    genome_index = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Sequence/bwa/UU_Cfam_GSD_1.0_ROSY.{amb,ann,bwt.2bit.64,pac,0123}"
    genome_dict = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Sequence/bwa/UU_Cfam_GSD_1.0_ROSY.dict"
    genome_fai = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Sequence/bwa/UU_Cfam_GSD_1.0_ROSY.fa.fai"
    known_variants = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Variation/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz"
    known_variants_index = "/data/references/Canis_lupus_familiaris/UCSC/UU_Cfam_GSD_1.0_canFam4/Variation/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz.tbi"

    // Tool parameters
    // Fastp parameters
    fastp_qualified_quality = 15        // More permissive for FFPE
    fastp_length_required = 50          // Min length for R1/R2
    fastp_umi_length = 8               // Expected UMI length
    
    // Boolean flags for trimming
    fastp_cut_front = true             // Enable 5' trimming for FFPE
    fastp_cut_tail = true              // Enable 3' trimming for FFPE
    fastp_cut_window_size = 4          // Window size for quality trimming
    fastp_cut_mean_quality = 20        // Quality threshold for trimming
    //aligner
    aligner = 'bwa-mem2'  // Options: 'bwa', 'bwa-mem2', 'dragen'
    min_read_length = 50
    min_base_quality = 20

    // UMI parameters
    umi_length = 8

    // Variant Calling parameters
    min_mapping_quality = 20

    // snpEff parameters
    snpeff_path = "/home/vjaganna/software/snpEff/exec/snpeff"
    snpeff_config = "/home/vjaganna/software/snpEff/snpEff.config"
    snpeff_memory = "4g"

    // Computing resources
    max_memory = '50.GB'
    max_cpus = 16
    max_time = '240.h'
}
log.info """
         UMI-FFPE VARIANT CALLING PIPELINE
         ================================
         genome: ${params.genome}
         reads : ${params.reads}
         outdir: ${params.outdir}
         aligner: ${params.aligner}
         """
workflow {
    // Input channels
    reads_ch = Channel
        .fromFilePairs(params.reads, size: 12) { file -> 
            def matcher = file =~ /(S[0-9]+Spleen)_L00[1-4]_R[1-3]_001\.fastq\.gz/
            matcher[0][1]
        }
        .map { sample_id, files -> 
            def grouped_files = files.groupBy { file ->
                def lane_matcher = file =~ /.*_(L00[1-4])_.*/
                lane_matcher[0][1]
            }
            def lanes = grouped_files.collectEntries { lane, lane_files ->
                [(lane): [
                    r1: lane_files.find { it.name.contains('_R1_') },
                    r2: lane_files.find { it.name.contains('_R2_') },
                    r3: lane_files.find { it.name.contains('_R3_') }
                ]]
            }
            [sample_id, lanes]
        }
        // Reference channels
        genome_ch = Channel.fromPath(params.genome)
        known_variants_ch = Channel.fromPath(params.known_variants)

        // Workflow execution
    TRIM_READS(reads_ch)
    ALIGN_READS(TRIM_READS.out.trimmed_reads, genome_ch)
    ANNOTATE_UMIS(ALIGN_READS.out.bam)
    MERGE_BAMS(ANNOTATE_UMIS.out.bam.groupTuple())
    MARK_DUPLICATES(MERGE_BAMS.out.merged_bam)
    BQSR(MARK_DUPLICATES.out.deduped_bam, genome_ch, known_variants_ch)
    CALL_VARIANTS(BQSR.out.recal_bam, genome_ch)
    GENOTYPE_GVCFS(CALL_VARIANTS.out.gvcf, genome_ch)
    FILTER_VARIANTS(GENOTYPE_GVCFS.out.vcf))
    ANNOTATE_VARIANTS(FILTER_VARIANTS.out.filtered_variants)
}

process TRIM_READS {
    tag "${sample_id}-${lane}"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
    tuple val(sample_id), val(lanes)
    
    output:
    tuple val(sample_id), val(lane), path("*_R{1,2,3}_trimmed.fastq.gz"), emit: trimmed_reads
    path "*.{html,json}", emit: trim_reports
    path "*.log", emit: logs
    
    script:
    def cut_front = params.fastp_cut_front ? '--cut_front' : ''
    def cut_tail = params.fastp_cut_tail ? '--cut_tail' : ''
    
    """
    # Process each lane
    for lane in "${lanes.keySet()}"; do
        # Process R1 and R2 (genomic reads)
        fastp \
            --in1 ${lanes[lane].r1} \
            --in2 ${lanes[lane].r2} \
            --out1 ${sample_id}_\${lane}_R1_trimmed.fastq.gz \
            --out2 ${sample_id}_\${lane}_R2_trimmed.fastq.gz \
            --detect_adapter_for_pe \
            --qualified_quality_phred ${params.fastp_qualified_quality} \
            --length_required ${params.fastp_length_required} \
            ${cut_front} \
            ${cut_tail} \
            --cut_window_size ${params.fastp_cut_window_size} \
            --cut_mean_quality ${params.fastp_cut_mean_quality} \
            --html ${sample_id}_\${lane}_R1R2.html \
            --json ${sample_id}_\${lane}_R1R2.json \
            --thread ${task.cpus} \
            --report_title "${sample_id}_\${lane}_R1R2" \
            --overrepresentation_analysis \
            2> ${sample_id}_\${lane}_R1R2.log

        # Process R3 (UMI read) separately
        fastp \
            --in1 ${lanes[lane].r3} \
            --out1 ${sample_id}_\${lane}_R3_trimmed.fastq.gz \
            --length_required ${params.fastp_umi_length} \
            --qualified_quality_phred ${params.fastp_qualified_quality} \
            ${cut_front} \
            ${cut_tail} \
            --cut_window_size ${params.fastp_cut_window_size} \
            --cut_mean_quality ${params.fastp_cut_mean_quality} \
            --html ${sample_id}_\${lane}_R3.html \
            --json ${sample_id}_\${lane}_R3.json \
            --thread ${task.cpus} \
            --report_title "${sample_id}_\${lane}_R3" \
            2> ${sample_id}_\${lane}_R3.log
    done
    """
}

process ALIGN_READS {
    tag "${sample_id}-${lane}"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(trimmed_reads)
    path genome
    
    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    path "*.bam.bai", emit: bai
    path "*_alignment_metrics.txt", emit: metrics
    
    script:
    def rgID = "${sample_id}"
    def platform = "ILLUMINA"
    def libraryID = "${sample_id}_lib"
    def assembly = "UU_Cfam_GSD_1.0"
    
    if (params.aligner == 'dragen-os') {
        """
        dragen-os \
            -r ${genome} \
            --bam-input ${trimmed_reads[0]} \
            --num-threads ${task.cpus} \
            --interleaved \
            > ${sample_id}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${sample_id}.bam ${sample_id}.sam
        samtools index ${sample_id}.bam
        rm ${sample_id}.sam

        # Collect alignment metrics
        samtools stats ${sample_id}.bam > ${sample_id}_alignment_metrics.txt
        """
    } else if (params.aligner == 'bwa-mem2') {
        """
        bwa-mem2 mem \
            -K 100000000 \
            -Y \
            -R "@RG\\tID:${rgID}\\tPL:${platform}\\tPU:${platform}\\tSM:${sample_id}\\tLB:${libraryID}\\tDS:${assembly}\\tCN:UNIBE" \
            -t ${task.cpus} \
            ${genome} \
            ${trimmed_reads[0]} \
            ${trimmed_reads[1]} \
            > ${sample_id}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${sample_id}.bam ${sample_id}.sam
        samtools index ${sample_id}.bam
        rm ${sample_id}.sam

        # Collect alignment metrics
        samtools stats ${sample_id}.bam > ${sample_id}_alignment_metrics.txt
        """
    } else {
        """
        bwa mem \
            -K 100000000 \
            -Y \
            -R "@RG\\tID:${rgID}\\tPL:${platform}\\tPU:${platform}\\tSM:${sample_id}\\tLB:${libraryID}\\tDS:${assembly}\\tCN:UNIBE" \
            -t ${task.cpus} \
            ${genome} \
            ${trimmed_reads[0]} \
            ${trimmed_reads[1]} \
            > ${sample_id}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${sample_id}.bam ${sample_id}.sam
        samtools index ${sample_id}.bam
        rm ${sample_id}.sam

        # Collect alignment metrics
        samtools stats ${sample_id}.bam > ${sample_id}_alignment_metrics.txt
        """
    }
}

process ANNOTATE_UMIS {
    tag "${sample_id}"
    publishDir "${params.outdir}/umi_annotated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("*_umi.bam"), path("*_umi.bam.bai"), emit: bam
    
    script:
    """
    # Add UMIs to BAM file using fgbio
    fgbio AnnotateBamWithUmis \
        -i ${bam} \
        -f ${sample_id}_R3_trimmed.fastq.gz \
        -o ${sample_id}_umi.bam
    
    # Index UMI-annotated BAM
    samtools index ${sample_id}_umi.bam
    """
}
process MERGE_BAMS {
    tag "${sample_id}"
    publishDir "${params.outdir}/merged", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bams), path(bais)
    
    output:
    tuple val(sample_id), path("*_merged.bam"), path("*_merged.bam.bai"), emit: merged_bam
    
    script:
    """
    # Merge lane-specific BAMs
    picard MergeSamFiles \
        ${bams.collect{ "I=$it" }.join(' ')} \
        O=${sample_id}_merged.bam \
        MERGE_SEQUENCE_DICTIONARIES=true \
        USE_THREADING=true \
        CREATE_INDEX=true
    """
}
process MARK_DUPLICATES {
    tag "${sample_id}"
    publishDir "${params.outdir}/deduped", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("*_deduped.bam"), path("*_deduped.bam.bai"), emit: deduped_bam
    path "*_duplicate_metrics.txt", emit: metrics
    
    script:
    """
    # UMI-aware duplicate marking
    picard UmiAwareMarkDuplicatesWithMateCigar \
        I=${bam} \
        O=${sample_id}_deduped.bam \
        M=${sample_id}_duplicate_metrics.txt \
        UMI_METRICS=${sample_id}_umi_metrics.txt \
        CREATE_INDEX=true \
        UMI_TAG_NAME=RX \
        ASSUME_SORTED=true
    """
}
process BQSR {
    tag "${sample_id}"
    publishDir "${params.outdir}/bqsr", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome
    path known_variants
    
    output:
    tuple val(sample_id), path("*_recal.bam"), path("*_recal.bam.bai"), emit: recal_bam
    path "*_recal.table", emit: metrics
    
    script:
    """
    # Generate recalibration table
    gatk BaseRecalibrator \
        -R ${genome} \
        -I ${bam} \
        --known-sites ${known_variants} \
        -O ${sample_id}_recal.table

    # Apply BQSR
    gatk ApplyBQSR \
        -R ${genome} \
        -I ${bam} \
        --bqsr-recal-file ${sample_id}_recal.table \
        -O ${sample_id}_recal.bam
    """
}
process CALL_VARIANTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path genome
    
    output:
    tuple val(sample_id), path("*_variants.vcf.gz"), path("*_variants.vcf.gz.tbi"), emit: variants
    
    script:
    """
    # Call variants using HaplotypeCaller
    gatk HaplotypeCaller \
        -R ${genome} \
        -I ${bam} \
        -O ${sample_id}_variants.vcf.gz \
        --native-pair-hmm-threads ${task.cpus} \
        -stand-call-conf 30 \
        --annotation QualByDepth \
        --annotation FisherStrand \
        --annotation RMSMappingQuality \
        --annotation ReadPosRankSumTest \
        -ERC GVCF

    # Index VCF
    gatk IndexFeatureFile -I ${sample_id}_variants.vcf.gz
    """
}
process GENOTYPE_GVCFS {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants/vcf", mode: 'copy'
    
    input:
    tuple val(sample_id), path(gvcf), path(tbi)
    path genome
    
    output:
    tuple val(sample_id), path("*_genotyped.vcf.gz"), path("*_genotyped.vcf.gz.tbi"), emit: vcf
    
    script:
    """
    # Genotype GVCF to produce final VCF
    gatk GenotypeGVCFs \
        -R ${genome} \
        -V ${gvcf} \
        -O ${sample_id}_genotyped.vcf.gz \
        --max-alternate-alleles 3 \
        -stand-call-conf 30

    # Index VCF
    gatk IndexFeatureFile -I ${sample_id}_genotyped.vcf.gz
    """
}
process FILTER_VARIANTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants/filtered", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    
    output:
    tuple val(sample_id), path("*_filtered.vcf.gz"), path("*_filtered.vcf.gz.tbi"), emit: filtered_variants
    
    script:
    """
    # Apply variant filters
    gatk VariantFiltration \
        -V ${vcf} \
        -O ${sample_id}_filtered.vcf.gz \
        --filter-name "QD_filter" --filter-expression "QD < 2.0" \
        --filter-name "FS_filter" --filter-expression "FS > 60.0" \
        --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
        --filter-name "Coverage_filter" --filter-expression "DP < 10"

    # Index filtered VCF
    gatk IndexFeatureFile -I ${sample_id}_filtered.vcf.gz
    """
}

process ANNOTATE_VARIANTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/variants/annotated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    
    output:
    tuple val(sample_id), path("*_annotated.vcf"), emit: annotated_vcf
    path "*_snpEff_stats.csv", emit: stats
    
    script:
    """
    ${params.snpeff_path} \
        -Xmx${params.snpeff_memory} \
        -csvStats ${sample_id}_snpEff_stats.csv \
        -c ${params.snpeff_config} \
        -v UU_Cfam_GSD_1.0_ROSY \
        ${vcf} \
        > ${sample_id}_annotated.vcf
    """
}


