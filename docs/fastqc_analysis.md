# FastQC Analysis

## Project

Microbiome Bioinformatics Training

Learner: Dr Guy Roussel Takuissu Nguemto

## Input data

File: data/example_reads.fastq.gz

Number of reads: 4

Read length: 12 bp

Encoding identified by FastQC: Illumina 1.5

## FastQC analysis

FastQC version: 0.12.1

The FASTQ file was analysed using FastQC v0.12.1.

### Basic Statistics

- Total sequences: 4
- Sequence length: 12 bp
- Overall GC content: 47%
- Encoding: Illumina 1.5

### Per base sequence quality

Status: FAIL

The quality profile was approximately Q9 across the 12 positions.

The profile was horizontal because all reads were generated with the same quality string.

### Per sequence quality scores

Status: FAIL

The sequence quality distribution was concentrated around Q9.

The low quality score is related to the artificial quality string used in this training dataset and the Illumina 1.5 encoding identified by FastQC.

### Per base sequence content

Status: FAIL

All four nucleotides (A, C, G and T) were present, but their proportions varied substantially between positions.

This reflects the highly simplified and artificial nature of the training sequences.

### Per sequence GC content

The individual reads had different GC contents, ranging approximately from 33% to 58%.

Because only four reads were available, the observed GC distribution was too small and discrete for meaningful biological interpretation.

### Sequence Length Distribution

Status: PASS

All four reads had a length of 12 bp.

### Sequence Duplication Levels

Status: PASS

The four reads had distinct nucleotide sequences. No sequence duplication was observed.

### Overrepresented sequences

No individual sequence was observed as overrepresented.

Interpretation should be cautious because the dataset contains only four reads.

### Adapter Content

Status: WARN

FastQC reported:

"Can't analyse adapters as read length is too short."

The reads were only 12 bp long, which is insufficient for meaningful adapter-content analysis.

This WARN should not be interpreted as evidence that adapters are present.

## Overall interpretation

This dataset is an artificial training dataset containing only four 12-bp reads. Therefore, several FastQC modules are not suitable for biological interpretation.

The purpose of this analysis was to learn how to perform and interpret basic FASTQ quality control rather than to assess the quality of a real sequencing library.

Important concepts learned include FASTQ structure, Phred quality scores, FASTQ quality encoding, nucleotide composition, GC content, sequence duplication, overrepresented sequences, and adapter-content assessment.

## Reproducibility

FastQC command:

fastqc data/example_reads.fastq.gz -o results/fastqc

FastQC version:

0.12.1

Output files:

- results/fastqc/example_reads_fastqc.html
- results/fastqc/example_reads_fastqc.zip
