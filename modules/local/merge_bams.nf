process MERGE_BAMS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::picard=2.25.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.25.5--hdfd78af_0' :
        'biocontainers/picard:2.25.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(bams), path(bais)  // [sample_id, [bam_files], [bai_files]]

    output:
    tuple val(meta), path("*_merged.bam"), path("*_merged.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def input_bams = bams.collect{ "I=$it" }.join(' ')
    """
    picard \\
        -Xmx${avail_mem}g \\
        MergeSamFiles \\
        $input_bams \\
        O=${prefix}_merged.bam \\
        MERGE_SEQUENCE_DICTIONARIES=true \\
        USE_THREADING=true \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MergeSamFiles --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
