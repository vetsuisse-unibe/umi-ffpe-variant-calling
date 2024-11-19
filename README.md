# UMI-based FFPE Variant Calling Pipeline

A Nextflow pipeline for processing UMI-based FFPE samples for variant calling.

## Pipeline Overview

The pipeline performs the following steps:
1. Read trimming and quality control (fastp)
2. Read alignment (BWA-MEM2/BWA/DRAGEN)
3. UMI annotation (fgbio)
4. BAM merging and duplicate marking
5. Base quality score recalibration (GATK)
6. Variant calling (GATK HaplotypeCaller)
7. Variant filtering and annotation (SnpEff)

## Usage

```bash
nextflow run main.nf -profile unibe,singularity
```

## Requirements

- Nextflow (>=21.10.3)
- Singularity/Docker/Conda
- Reference genome and known variants
