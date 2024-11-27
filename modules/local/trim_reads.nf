process TRIM_READS {
    tag "$meta.id"

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path "*.{html,json}", emit: reports
    path "*.log", emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash

    # Group files by lane and read number
    declare -A r1_files r2_files r3_files
    
    for file in ${reads}; do
        if [[ \$file =~ L00([0-9])_R([0-9])_001.fastq.gz ]]; then
            lane="\${BASH_REMATCH[1]}"
            read="\${BASH_REMATCH[2]}"
            
            case \$read in
                1) r1_files[\$lane]=\$file ;;  # Template Read 1
                2) r2_files[\$lane]=\$file ;;  # UMI Read
                3) r3_files[\$lane]=\$file ;;  # Template Read 2
            esac
        fi
    done

    # Process each lane
    for lane in \${!r1_files[@]}; do
        echo "Processing lane \$lane:"
        echo "Template R1: \${r1_files[\$lane]}"
        echo "UMI R2: \${r2_files[\$lane]}"
        echo "Template R3: \${r3_files[\$lane]}"

        # Process Template reads (R1 and R3)
        fastp \\
            --in1 \${r1_files[\$lane]} \\
            --in2 \${r3_files[\$lane]} \\
            --out1 ${prefix}_L00\${lane}_R1_trimmed.fastq.gz \\
            --out2 ${prefix}_L00\${lane}_R3_trimmed.fastq.gz \\
            --detect_adapter_for_pe \\
            --qualified_quality_phred 15 \\
            --length_required 50 \\
            --cut_front \\
            --cut_tail \\
            --cut_window_size 4 \\
            --cut_mean_quality 15 \\
            --html ${prefix}_L00\${lane}_template.html \\
            --json ${prefix}_L00\${lane}_template.json \\
            --thread $task.cpus \\
            --report_title "${prefix}_L00\${lane}_template" \\
            --overrepresentation_analysis \\
            2> ${prefix}_L00\${lane}_template.log

        # Process UMI read (R2) separately
        fastp \\
            --in1 \${r2_files[\$lane]} \\
            --out1 ${prefix}_L00\${lane}_R2_trimmed.fastq.gz \\
            --length_required 8 \\
            --qualified_quality_phred 15 \\
            --cut_front \\
            --cut_tail \\
            --cut_window_size 4 \\
            --cut_mean_quality 15 \\
            --html ${prefix}_L00\${lane}_umi.html \\
            --json ${prefix}_L00\${lane}_umi.json \\
            --thread $task.cpus \\
            --report_title "${prefix}_L00\${lane}_umi" \\
            2> ${prefix}_L00\${lane}_umi.log
    done

    cat <<-END_VERSIONS > versions.yml
    "TRIM_READS":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}