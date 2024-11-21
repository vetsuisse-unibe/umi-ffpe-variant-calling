#!/bin/bash

#SBATCH --mail-user=vidhya.jagannathan@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name=variant_calling
#SBATCH --output=pipeline_run_%j.log
#SBATCH --error=pipeline_run_%j.err
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --partition=pibu_el8


# Load the required modules
module load Java/17.0.6

# Run the pipeline
nextflow run workflows/main.nf \
    -profile unibe,singularity \
    -resume \
    -with-report execution_report.html \
    -with-trace execution_trace.txt \
    -with-timeline timeline.html \
    -with-dag flowchart.png