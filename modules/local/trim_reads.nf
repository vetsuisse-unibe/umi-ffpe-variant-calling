// modules/local/trim_reads.nf

process TRIM_READS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), val(lanes)  // [sample_id, [lane: [r1:file, r2:file, r3:file]]]

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path "*.{html,json}", emit: reports
    path "*.log", emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Process each lane
    for lane in "${lanes.keySet()}"; do
        # Process R1 and R2 (genomic reads)
        fastp \\
            --in1 ${lanes[lane].r1} \\
            --in2 ${lanes[lane].r2} \\
            --out1 ${prefix}_\${lane}_R1_trimmed.fastq.gz \\
            --out2 ${prefix}_\${lane}_R2_trimmed.fastq.gz \\
            --detect_adapter_for_pe \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --length_required ${params.fastp_length_required} \\
            ${params.fastp_cut_front ? '--cut_front' : ''} \\
            ${params.fastp_cut_tail ? '--cut_tail' : ''} \\
            --cut_window_size ${params.fastp_cut_window_size} \\
            --cut_mean_quality ${params.fastp_cut_mean_quality} \\
            --html ${prefix}_\${lane}_R1R2.html \\
            --json ${prefix}_\${lane}_R1R2.json \\
            --thread $task.cpus \\
            --report_title "${prefix}_\${lane}_R1R2" \\
            --overrepresentation_analysis \\
            $args \\
            2> ${prefix}_\${lane}_R1R2.log

        # Process R3 (UMI read) separately
        fastp \\
            --in1 ${lanes[lane].r3} \\
            --out1 ${prefix}_\${lane}_R3_trimmed.fastq.gz \\
            --length_required ${params.fastp_umi_length} \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            ${params.fastp_cut_front ? '--cut_front' : ''} \\
            ${params.fastp_cut_tail ? '--cut_tail' : ''} \\
            --cut_window_size ${params.fastp_cut_window_size} \\
            --cut_mean_quality ${params.fastp_cut_mean_quality} \\
            --html ${prefix}_\${lane}_R3.html \\
            --json ${prefix}_\${lane}_R3.json \\
            --thread $task.cpus \\
            --report_title "${prefix}_\${lane}_R3" \\
            $args \\
            2> ${prefix}_\${lane}_R3.log
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}