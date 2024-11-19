process MARK_DUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=2.25.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.5--hdfd78af_0' :
        'biocontainers/picard:2.25.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)  // [sample_id, bam, bai]

    output:
    tuple val(meta), path("*_deduped.bam"), path("*_deduped.bam.bai"), emit: bam
    path "*_duplicate_metrics.txt", emit: metrics
    path "*_umi_metrics.txt", emit: umi_metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        -Xmx${avail_mem}g \\
        UmiAwareMarkDuplicatesWithMateCigar \\
        I=${bam} \\
        O=${prefix}_deduped.bam \\
        M=${prefix}_duplicate_metrics.txt \\
        UMI_METRICS=${prefix}_umi_metrics.txt \\
        CREATE_INDEX=true \\
        UMI_TAG_NAME=RX \\
        ASSUME_SORTED=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard UmiAwareMarkDuplicatesWithMateCigar --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
