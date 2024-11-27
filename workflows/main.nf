#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIM_READS } from '../modules/local/trim_reads'
include { ALIGN_READS } from '../modules/local/align_reads'

workflow {
    // Input channels
    Channel
        .fromFilePairs("/data/projects/p525_dog_wgs/nextflow/ffpe/raw_data/fastq/*_L00{1,2,3,4}_R{1,2,3}_001.fastq.gz", size: -1)
        .map { id, files -> 
            def meta = [id: id]
            return [meta, files]
        }
        .set { ch_reads }

    // Step 1: Trim reads
    TRIM_READS(ch_reads)
}
