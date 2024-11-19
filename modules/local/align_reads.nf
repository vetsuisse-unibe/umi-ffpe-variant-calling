// modules/local/align_reads.nf

process ALIGN_READS {
    tag "$meta.id"
    label 'process_high'

    conda (params.aligner == 'bwa-mem2') ? "bioconda::bwa-mem2=2.2.1" : "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        (params.aligner == 'bwa-mem2') ? 'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he513fc3_0' :
        (params.aligner == 'bwa') ? 'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'https://depot.galaxyproject.org/singularity/dragmap:1.3.0--h91baf5a_3' :
        (params.aligner == 'bwa-mem2') ? 'biocontainers/bwa-mem2:2.2.1--he513fc3_0' :
        (params.aligner == 'bwa') ? 'biocontainers/bwa:0.7.17--hed695b0_7' :
        'biocontainers/dragmap:1.3.0--h91baf5a_3' }"

    input:
    tuple val(meta), path(reads)  // [sample_id, [trimmed_fastq_files]]
    path genome                   // Reference genome fasta

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam
    path "*_alignment_metrics.txt", emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rgID = prefix
    def platform = "ILLUMINA"
    def libraryID = "${prefix}_lib"
    def assembly = "UU_Cfam_GSD_1.0"
    
    if (params.aligner == 'dragen-os') {
        """
        dragen-os \\
            -r ${genome} \\
            --bam-input ${reads[0]} \\
            --num-threads ${task.cpus} \\
            --interleaved \\
            $args \\
            > ${prefix}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${prefix}.bam ${prefix}.sam
        samtools index ${prefix}.bam
        rm ${prefix}.sam

        samtools stats ${prefix}.bam > ${prefix}_alignment_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dragen-os: \$(dragen-os --version 2>&1)
            samtools: \$(samtools --version | grep ^samtools | sed 's/^.*samtools //; s/ .*\$//')
        END_VERSIONS
        """
    } else if (params.aligner == 'bwa-mem2') {
        """
        bwa-mem2 mem \\
            -K 100000000 \\
            -Y \\
            -R "@RG\\tID:${rgID}\\tPL:${platform}\\tPU:${platform}\\tSM:${prefix}\\tLB:${libraryID}\\tDS:${assembly}\\tCN:UNIBE" \\
            -t ${task.cpus} \\
            $args \\
            ${genome} \\
            ${reads[0]} \\
            ${reads[1]} \\
            > ${prefix}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${prefix}.bam ${prefix}.sam
        samtools index ${prefix}.bam
        rm ${prefix}.sam

        samtools stats ${prefix}.bam > ${prefix}_alignment_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa-mem2: \$(bwa-mem2 version 2>&1 | sed 's/^.*bwa-mem2 //; s/Using.*\$//')
            samtools: \$(samtools --version | grep ^samtools | sed 's/^.*samtools //; s/ .*\$//')
        END_VERSIONS
        """
    } else {
        """
        bwa mem \\
            -K 100000000 \\
            -Y \\
            -R "@RG\\tID:${rgID}\\tPL:${platform}\\tPU:${platform}\\tSM:${prefix}\\tLB:${libraryID}\\tDS:${assembly}\\tCN:UNIBE" \\
            -t ${task.cpus} \\
            $args \\
            ${genome} \\
            ${reads[0]} \\
            ${reads[1]} \\
            > ${prefix}.sam

        samtools sort -m 5G -@ ${task.cpus} -o ${prefix}.bam ${prefix}.sam
        samtools index ${prefix}.bam
        rm ${prefix}.sam

        samtools stats ${prefix}.bam > ${prefix}_alignment_metrics.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwa: \$(bwa 2>&1 | grep ^Version | sed 's/^.*Version: //; s/$//')
            samtools: \$(samtools --version | grep ^samtools | sed 's/^.*samtools //; s/ .*\$//')
        END_VERSIONS
        """
    }
}