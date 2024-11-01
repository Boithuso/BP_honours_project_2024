#!/bin/bash

#the input directory must contain bam files that are indexed and have there read groups
#within the gatk directory my have a list "input_bams.list" with the bam files

python3 gatk DepthOfCoverage -R genome.fa -O hla_DOC -I input_bams.list \
-DF NotSecondaryAlignmentReadFilter \
-gene-list hla_genes_and_gene_predictions.refseq \
-L GL000250.2 -L GL000251.2 -L GL000252.2 -L GL000253.2 \
-L GL000254.2 -L GL000255.2 -L GL000256.2 \
-L chr6:29700000-33500000
