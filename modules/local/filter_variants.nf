process FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.2.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0' :
        'biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // [sample_id, vcf, index]

    output:
    tuple val(meta), path("*_filtered.vcf.gz"), path("*_filtered.vcf.gz.tbi"), emit: vcf
    path "*_filtered.vcf.gz.stats", optional: true, emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK VariantFiltration] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}g" VariantFiltration \\
        -V ${vcf} \\
        -O ${prefix}_filtered.vcf.gz \\
        --tmp-dir . \\
        --filter-name "QD_filter" --filter-expression "QD < ${params.qd_threshold ?: '2.0'}" \\
        --filter-name "FS_filter" --filter-expression "FS > ${params.fs_threshold ?: '60.0'}" \\
        --filter-name "MQ_filter" --filter-expression "MQ < ${params.mq_threshold ?: '40.0'}" \\
        --filter-name "SOR_filter" --filter-expression "SOR > ${params.sor_threshold ?: '3.0'}" \\
        --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < ${params.rprs_threshold ?: '-8.0'}" \\
        --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < ${params.mqrs_threshold ?: '-12.5'}" \\
        --filter-name "Coverage_filter" --filter-expression "DP < ${params.min_coverage ?: '10'}" \\
        $args

    # Generate filter stats if requested
    if [ "${params.generate_filter_stats}" = "true" ]; then
        gatk VariantEval \\
            -V ${prefix}_filtered.vcf.gz \\
            -O ${prefix}_filtered.vcf.gz.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
