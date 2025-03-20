#!/bin/bash

mkdir -p DNA_transposon_fasta
mkdir -p genomes
node tableify.ts | while read bed twobit accession; do
  echo $bed $twobit $accession

  # Store basename in a variable
  twobit_base=$(basename $twobit)
  bed_base=$(basename $bed .bed)

  wget -nc -q -P genomes $twobit
  twoBitToFa genomes/$twobit_base genomes/$twobit_base.fa
  samtools faidx genomes/$twobit_base.fa

  bedtools slop -i DNA_transposons/$bed -g genomes/$twobit_base.fa.fai -b 2000 |
    gzip >DNA_transposons/$bed_base.slop.bed.gz

  ## output slop'd fa
  bedtools getfasta -bed DNA_transposons/$bed_base.slop.bed.gz -fi genomes/$twobit_base.fa -name+ |
    bgzip >DNA_transposon_fasta/$accession.dna_transposons.slop.fa.gz

  ## output non-slop'd fa
  bedtools getfasta -bed DNA_transposons/$bed -fi genomes/$twobit_base.fa -name+ |
    bgzip >DNA_transposon_fasta/$accession.dna_transposons.fa.gz

  rm -f genomes/$twobit_base.fa*

done
