// modules/local/bqsr.nf

process BQSR {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.2.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0' :
        'biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)     // [sample_id, bam, bai]
    path(genome)                               // Reference genome fasta
    path(known_sites)                          // Known variants for recalibration

    output:
    tuple val(meta), path("*_recal.bam"), path("*_recal.bam.bai"), emit: bam
    path "*_recal.table", emit: table
    path "*_pre_recal.table", optional: true, emit: pre_table
    path "*_recalibration_plots.pdf", optional: true, emit: recal_plots
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def known_sites_command = known_sites ? "--known-sites ${known_sites}" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    def analyze_covariates_command = params.generate_plots ? """
        gatk AnalyzeCovariates \\
            -before ${prefix}_pre_recal.table \\
            -after ${prefix}_recal.table \\
            -plots ${prefix}_recalibration_plots.pdf
    """ : ""

    """
    # First pass of recalibration
    gatk --java-options "-Xmx${avail_mem}g" BaseRecalibrator \\
        -R ${genome} \\
        -I ${bam} \\
        ${known_sites_command} \\
        -O ${prefix}_recal.table \\
        --tmp-dir . \\
        $args

    # Apply BQSR
    gatk --java-options "-Xmx${avail_mem}g" ApplyBQSR \\
        -R ${genome} \\
        -I ${bam} \\
        --bqsr-recal-file ${prefix}_recal.table \\
        -O ${prefix}_recal.bam \\
        --tmp-dir . \\
        $args

    if [ "${params.validate_bqsr}" == "true" ]; then
        # Second pass of recalibration for validation
        gatk --java-options "-Xmx${avail_mem}g" BaseRecalibrator \\
            -R ${genome} \\
            -I ${prefix}_recal.bam \\
            ${known_sites_command} \\
            -O ${prefix}_post_recal.table \\
            --tmp-dir . \\
            $args

        # Generate recalibration plots
        $analyze_covariates_command
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}