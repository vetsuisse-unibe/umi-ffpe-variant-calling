process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::snpeff=5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_0' :
        'biocontainers/snpeff:5.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_index)  // [sample_id, vcf, index]

    output:
    tuple val(meta), path("*_annotated.vcf.gz"), path("*_annotated.vcf.gz.tbi"), emit: vcf
    path "*_snpEff_stats.csv", emit: stats
    path "*_snpEff.csv", emit: report
    path "*_snpEff.genes.txt", emit: genes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 3
    if (!task.memory) {
        log.info '[SnpEff] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def snpeff_db = params.snpeff_db ?: 'UU_Cfam_GSD_1.0_ROSY'
    
    """
    # Run SnpEff
    snpEff \\
        -Xmx${avail_mem}g \\
        -v \\
        -stats ${prefix}_snpEff_stats.csv \\
        -csvStats ${prefix}_snpEff.csv \\
        -noLog \\
        ${args} \\
        ${snpeff_db} \\
        ${vcf} > ${prefix}_annotated.vcf

    # Compress and index output VCF
    bgzip -f ${prefix}_annotated.vcf
    tabix -p vcf ${prefix}_annotated.vcf.gz

    # Extract gene-level summary
    cat ${prefix}_snpEff.csv | grep -v '^#' | cut -f1,2,3 | sort | uniq > ${prefix}_snpEff.genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -d ' ' -f2)
        bgzip: \$(echo \$(bgzip -h 2>&1) | grep Version | sed 's/^.*Version: //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | grep Version | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
