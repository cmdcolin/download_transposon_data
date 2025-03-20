#!/bin/bash
rm -rf extracts
mkdir -p extracts

for i in DNA_transposon_fasta/*.fa.gz; do
  echo $i
  IFS='_' read -ra parts2 <<<"$i"
  X=$(basename ${parts2[1]} .dna)
  cat $i.fai | while read r x y z; do
    IFS='::' read -a parts <<<"$r"

    samtools faidx $i $r --mark-strand custom," ${parts2[0]}_$X" >>extracts/$parts.fa
  done
done
