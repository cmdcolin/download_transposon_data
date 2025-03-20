#!/bin/bash

mkdir DNA_transposon_fasta
node tableify.ts | while read bed twobit accession; do
  echo $bed $twobit $accession
  twoBitToFa $twobit $(basename $twobit).fa
  samtools faidx $(basename $twobit).fa

  bedtools slop -i DNA_transposons/$bed -g $(basename $twobit).fa.fai -l 2000 -r 2000 | gzip >DNA_transposon_fasta/slop.$bed.gz
  bedtools getfasta -bed DNA_transposon_fasta/slop.$bed.gz -fi $(basename $twobit).fa -name+ | bgzip >DNA_transposon_fasta/$accession.dna_transposons.fa.gz

  rm -f $(basename $twobit).fa*
done
