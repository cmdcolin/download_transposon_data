node tableify.ts | head -n1 | while read bed twobit accession; do
  echo $bed $twobit $accession
  wget -nc $twobit
  twoBitToFa $twobit $(basename $twobit).fa

  bedtools getfasta -bed DNA_transposons/$bed -fi $(basename $twobit).fa -name+ | bgzip >DNA_transposon_fastas/$accession.dnate.fa.gz
  rm -f $twobit $(basename $twobit).fa
done
