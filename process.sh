#!/bin/bash

mkdir -p DNA_transposon_fasta
node tableify.ts | while read bed twobit accession; do
  echo $bed $twobit $accession
  wget -nc -q -P genomes $twobit
  twoBitToFa genomes/$(basename twobit) genomes/$(basename $twobit).fa
  samtools faidx genomes/$(basename $twobit).fa

  bedtools slop -i DNA_transposons/$bed -g genomes/$(basename $twobit).fa.fai -b 2000 | gzip >DNA_transposons/slop.$bed.gz

  ## output slop'd fa
  bedtools getfasta -bed DNA_transposons/slop.$bed.gz -fi genomes/$(basename $twobit).fa -name+ | bgzip >DNA_transposon_fasta/$accession.dna_transposons.slop.fa.gz

  ## output non-slop'd fa
  bedtools getfasta -bed DNA_transposons/$bed -fi $(basename $twobit).fa -name+ | bgzip >DNA_transposon_fasta/$accession.dna_transposons.fa.gz

done
