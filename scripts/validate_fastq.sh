#!/bin/bash

FASTQ="$1"

if [ -z "$FASTQ" ]; then
    echo "Usage: $0 <fastq_file>"
    exit 1
fi

if [ ! -f "$FASTQ" ]; then
    echo "ERROR: File not found: $FASTQ"
    exit 1
fi

echo "======================================"
echo "FASTQ QC SUMMARY"
echo "======================================"
echo "File: $FASTQ"

LINES=$(wc -l < "$FASTQ")

echo "Total lines: $LINES"

if [ $((LINES % 4)) -ne 0 ]; then
    echo "ERROR: Number of lines is not divisible by 4."
    echo "STATUS: FAIL"
    exit 1
fi

READS=$((LINES / 4))
echo "Reads: $READS"

LENGTH_ERRORS=$(awk '
NR % 4 == 1 {id=$0}
NR % 4 == 2 {seqlen=length($0)}
NR % 4 == 0 {
    quallen=length($0)
    if (seqlen != quallen) count++
}
END {print count+0}
' "$FASTQ")

INVALID_SEQUENCES=$(awk '
NR % 4 == 2 && $0 !~ /^[ACGTN]+$/ {count++}
END {print count+0}
' "$FASTQ")

echo "Invalid lengths: $LENGTH_ERRORS"
echo "Invalid sequences: $INVALID_SEQUENCES"

echo "--------------------------------------"

if [ "$LENGTH_ERRORS" -eq 0 ] && [ "$INVALID_SEQUENCES" -eq 0 ]; then
    echo "STATUS: PASS"
    exit 0
else
    echo "STATUS: FAIL"
    exit 1
fi
