// modules/local/call_variants.nf

process CALL_VARIANTS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.2.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0' :
        'biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)     // [sample_id, bam, bai]
    path(genome)                               // Reference genome fasta

    output:
    tuple val(meta), path("*_variants.g.vcf.gz"), path("*_variants.g.vcf.gz.tbi"), emit: gvcf
    tuple val(meta), path("*_variants.vcf.gz"), path("*_variants.vcf.gz.tbi"), optional: true, emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = params.intervals ? "--intervals ${params.intervals}" : ""
    def emit_ref_confidence = params.emit_ref_confidence ? "--emit-ref-confidence ${params.emit_ref_confidence}" : "--emit-ref-confidence GVCF"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}g" HaplotypeCaller \\
        -R ${genome} \\
        -I ${bam} \\
        -O ${prefix}_variants.g.vcf.gz \\
        ${interval_command} \\
        ${emit_ref_confidence} \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir . \\
        -stand-call-conf ${params.stand_call_conf ?: '30'} \\
        --annotation QualByDepth \\
        --annotation FisherStrand \\
        --annotation RMSMappingQuality \\
        --annotation ReadPosRankSumTest \\
        --annotation MappingQualityRankSumTest \\
        --annotation StrandOddsRatio \\
        --annotation DepthPerAlleleBySample \\
        $args

    # Index the GVCF
    gatk IndexFeatureFile -I ${prefix}_variants.g.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}