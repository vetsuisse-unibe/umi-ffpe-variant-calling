#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//variant calling sub-workflow
include { VARIANT_CALLING } from './sub-workflows/variant_calling'


//Parameters
params {
    // Input/Output
    reads = "/data/projects/p525_dog_wgs/nextflow/ffpe/raw_data/fastq/*_L00{1,2,3,4}_R{1,2,3}_001.fastq.gz"
    outdir = "/data/projects/p525_dog_wgs/nextflow/ffpe/results"

    //Reference files
    genome = null
    genome_index = null
    genome_dict = null
    genome_fai = null
    known_variants = null
    known_variants_index = null

    // Tool parameters
    fastp_qualified_quality = 15
    fastp_length_required = 50
    fastp_umi_length = 8
    fastp_cut_front = true
    fastp_cut_tail = true
    fastp_cut_window_size = 4
    fastp_cut_mean_quality = 20
    aligner = 'bwa-mem2'
    min_read_length = 50
    min_base_quality = 20
    umi_length = 8
    min_mapping_quality = 20
    snpeff_path = null
    snpeff_config = null
    snpeff_memory = "4g"

    // Computing resources
    max_memory = '50.GB'
    max_cpus = 16
    max_time = '240.h'
}
// Load base.config by default for all pipelines
includeConfig 'conf/nextflow.config'

// Load modules.config for module specific options
includeConfig 'conf/modules.config'

// Load profiles.config
includeConfig 'conf/profiles.config'

// Function to validate required parameters
def validateParameters() {
    if (!params.genome) {
        exit 1, "Reference genome not specified!"
    }
    // Add other parameter validations as needed
}

// Main workflow

workflow {
    // Validate parameters
    validateParameters()

   // Create input channel
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

         // Run variant calling workflow
    VARIANT_CALLING(
        reads_ch,
        genome_ch,
        known_variants_ch
    )

}
 
 // Workflow completion notification
workflow.onComplete {
    log.info """
    Pipeline completed at: ${workflow.complete}
    Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Execution duration: ${workflow.duration}
    """
}
// Manifest
manifest {
    name = 'umi-ffpe-variant-calling'
    author = 'Vidhya Jagannathan'
    homePage = 'https://github.com/vetsuisse-unibe/umi-ffpe-variant-calling.git'
    description = 'Nextflow pipeline for UMI-based FFPE variant calling'
    mainScript = 'main.nf'
    nextflowVersion = '>=21.10.3'
    version = '1.0.0'
}
