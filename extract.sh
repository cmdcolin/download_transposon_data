#!/bin/bash
rm -rf extracts
mkdir -p extracts

for i in *.fa.gz; do
  echo $i
  cat $i.fai | while read r x y z; do
    IFS='::' read -a parts <<<"$r"
    IFS='_' read -ra parts2 <<<"$i"
    X=$(basename ${parts2[1]} .dna)
    samtools faidx $i $r --mark-strand custom," ${parts2[0]}_$X" >>extracts/$parts.fa
  done
done
