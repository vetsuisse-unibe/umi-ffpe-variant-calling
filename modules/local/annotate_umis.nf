// modules/local/annotate_umis.nf

process ANNOTATE_UMIS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fgbio=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.3.0--hdfd78af_0' :
        'biocontainers/fgbio:2.3.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)     // [sample_id, aligned_bam, bam_index]

    output:
    tuple val(meta), path("*_umi.bam"), path("*_umi.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[fgbio AnnotateBamWithUmis] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    fgbio \\
        --tmp-dir . \\
        AnnotateBamWithUmis \\
        -i ${bam} \\
        -f ${prefix}_R3_trimmed.fastq.gz \\
        -o ${prefix}_umi.bam \\
        $args

    samtools index ${prefix}_umi.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$(fgbio --version 2>&1 | sed 's/^Version: //')
        samtools: \$(samtools --version | grep ^samtools | sed 's/^.*samtools //; s/ .*\$//')
    END_VERSIONS
    """
}