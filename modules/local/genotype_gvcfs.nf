// modules/local/genotype_gvcfs.nf

process GENOTYPE_GVCFS {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::gatk4=4.2.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.1--hdfd78af_0' :
        'biocontainers/gatk4:4.2.6.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(gvcf), path(gvcf_index)  // [sample_id, gvcf, index]
    path(genome)                                    // Reference genome fasta

    output:
    tuple val(meta), path("*_genotyped.vcf.gz"), path("*_genotyped.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*_genotyped.vcf.gz.stats"), optional: true, emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def interval_command = params.intervals ? "--intervals ${params.intervals}" : ""
    def dbsnp_command = params.dbsnp ? "--dbsnp ${params.dbsnp}" : ""
    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK GenotypeGVCFs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    
    """
    gatk --java-options "-Xmx${avail_mem}g" GenotypeGVCFs \\
        -R ${genome} \\
        -V ${gvcf} \\
        -O ${prefix}_genotyped.vcf.gz \\
        ${interval_command} \\
        ${dbsnp_command} \\
        --tmp-dir . \\
        --max-alternate-alleles ${params.max_alternate_alleles ?: '3'} \\
        -stand-call-conf ${params.stand_call_conf ?: '30'} \\
        --annotation QualByDepth \\
        --annotation FisherStrand \\
        --annotation RMSMappingQuality \\
        --annotation ReadPosRankSumTest \\
        --annotation MappingQualityRankSumTest \\
        --annotation StrandOddsRatio \\
        --annotation DepthPerAlleleBySample \\
        $args

    # Generate VCF stats if requested
    if [ "${params.generate_vcf_stats}" = "true" ]; then
        gatk VariantEval \\
            -R ${genome} \\
            --eval ${prefix}_genotyped.vcf.gz \\
            ${dbsnp_command} \\
            -O ${prefix}_genotyped.vcf.gz.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}