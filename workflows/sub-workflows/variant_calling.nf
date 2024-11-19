#!/usr/bin/env nextflow

include { TRIM_READS } from '../../modules/local/trim_reads'
include { ALIGN_READS } from '../../modules/local/align_reads'
include { ANNOTATE_UMIS } from '../../modules/local/annotate_umis'
include { MERGE_BAMS } from '../../modules/local/merge_bams'
include { MARK_DUPLICATES } from '../../modules/local/mark_duplicates'
include { BQSR } from '../../modules/local/bqsr'
include { CALL_VARIANTS } from '../../modules/local/call_variants'
include { GENOTYPE_GVCFS } from '../../modules/local/genotype_gvcfs'
include { FILTER_VARIANTS } from '../../modules/local/filter_variants'
include { ANNOTATE_VARIANTS } from '../../modules/local/annotate_variants'

workflow VARIANT_CALLING {
    take:
    reads_ch       // channel: [ val(meta), [ reads ] ]
    genome_ch      // channel: /path/to/genome.fa
    known_vars_ch  // channel: /path/to/known_variants.vcf

    main:
    versions_ch = Channel.empty()

    // Trim reads
    TRIM_READS ( reads_ch )
    versions_ch = versions_ch.mix(TRIM_READS.out.versions)

    // Align reads
    ALIGN_READS (
        TRIM_READS.out.reads,
        genome_ch
    )
    versions_ch = versions_ch.mix(ALIGN_READS.out.versions)

    // Annotate UMIs
    ANNOTATE_UMIS ( ALIGN_READS.out.bam )
    versions_ch = versions_ch.mix(ANNOTATE_UMIS.out.versions)

    // Merge BAMs from multiple lanes
    MERGE_BAMS ( ANNOTATE_UMIS.out.bam.groupTuple() )
    versions_ch = versions_ch.mix(MERGE_BAMS.out.versions)

    // Mark duplicates
    MARK_DUPLICATES ( MERGE_BAMS.out.bam )
    versions_ch = versions_ch.mix(MARK_DUPLICATES.out.versions)

    // Base quality score recalibration
    BQSR (
        MARK_DUPLICATES.out.bam,
        genome_ch,
        known_vars_ch
    )
    versions_ch = versions_ch.mix(BQSR.out.versions)

    // Variant calling
    CALL_VARIANTS (
        BQSR.out.bam,
        genome_ch
    )
    versions_ch = versions_ch.mix(CALL_VARIANTS.out.versions)

    // Genotype GVCFs
    GENOTYPE_GVCFS (
        CALL_VARIANTS.out.gvcf,
        genome_ch
    )
    versions_ch = versions_ch.mix(GENOTYPE_GVCFS.out.versions)

    // Filter variants
    FILTER_VARIANTS ( GENOTYPE_GVCFS.out.vcf )
    versions_ch = versions_ch.mix(FILTER_VARIANTS.out.versions)

    // Annotate variants
    ANNOTATE_VARIANTS ( FILTER_VARIANTS.out.vcf )
    versions_ch = versions_ch.mix(ANNOTATE_VARIANTS.out.versions)

    emit:
    reads = TRIM_READS.out.reads                  // channel: [ val(meta), [ reads ] ]
    bam = ALIGN_READS.out.bam                     // channel: [ val(meta), [ bam ] ]
    recal_bam = BQSR.out.bam                      // channel: [ val(meta), [ bam ] ]
    variants = ANNOTATE_VARIANTS.out.vcf          // channel: [ val(meta), [ vcf ] ]
    versions = versions_ch                         // channel: [ versions.yml ]
}