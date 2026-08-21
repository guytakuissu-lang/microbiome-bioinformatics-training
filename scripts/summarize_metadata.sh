#!/bin/bash

echo "Total number of lines:"
wc -l data/example_metadata.tsv

echo "T2D samples:"
grep -w T2D data/example_metadata.tsv

echo "Sample identifiers and groups:"
cut -f1,2 data/example_metadata.tsv
