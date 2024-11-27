#!/bin/bash

#SBATCH --mail-user=vidhya.jagannathan@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name=trim_reads
#SBATCH --output=trim_reads_%j.log
#SBATCH --error=trim_reads_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=pibu_el8
#SBATCH --account=p525

# Load the required modules
module load Java/17.0.6

# Run the pipeline with only trim_reads
nextflow run workflows/main.nf \
    -profile unibe,singularity \
    -resume \
    -with-report trim_reads_report.html \
    -with-trace trim_reads_trace.txt \
    -with-timeline trim_reads_timeline.html \
    -with-dag trim_reads_flowchart.png