# MultiQC Analysis

## Project

Microbiome Bioinformatics Training

Learner: Dr Guy Roussel Takuissu Nguemto

## Objective

The objective was to aggregate and summarize FastQC quality-control results using MultiQC.

## Input data

Input FASTQ file:

data/example_reads.fastq.gz

Number of reads:

4

Read length:

12 bp

## FastQC results

FastQC version:

0.12.1

The FASTQ file was previously analysed using FastQC.

The main observations were:

- Total sequences: 4
- Sequence length: 12 bp
- Overall GC content: 47%
- Per base sequence quality: FAIL
- Per sequence quality scores: FAIL
- Per base sequence content: FAIL
- Sequence length distribution: PASS
- Sequence duplication levels: PASS
- Overrepresented sequences: no overrepresented sequence detected
- Adapter content: WARN because the reads were too short for meaningful adapter analysis

## MultiQC analysis

MultiQC version:

1.35

MultiQC was used to aggregate the FastQC results into a single HTML report.

Command used:

multiqc results/fastqc/ -o results/multiqc

Output directory:

results/multiqc/

Main output:

results/multiqc/multiqc_report.html

## Interpretation

The MultiQC report successfully summarized the FastQC results.

Because the training dataset contains only four artificial reads of 12 bp, the quality-control results should not be interpreted as representative of a real sequencing library.

The FAIL status for several FastQC modules is mainly related to the artificial quality scores, nucleotide composition and extremely small dataset size.

The WARN status for Adapter Content does not demonstrate the presence of adapters. FastQC could not perform a meaningful adapter-content analysis because the reads were only 12 bp long.

The main purpose of this exercise was to learn how MultiQC can aggregate and summarize sequencing quality-control results in a reproducible workflow.

## Reproducibility

FastQC command:

fastqc data/example_reads.fastq.gz -o results/fastqc

MultiQC command:

multiqc results/fastqc/ -o results/multiqc

FastQC version:

0.12.1

MultiQC version:

1.35
