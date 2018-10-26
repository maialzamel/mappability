#!/bin/bash
mkdir output
mkdir index
mkdir outputdump
k=30
e=1
o=10

for ((i=0; i<=31; i++)); do
./create_index -G ./data/$i.fasta -I ./index/$i
./mappability -I ./index/$i -K $k -E $e -o $o -O ./output/$i
./mappability_dump -I ./output/$i'_'$e'_'$k'_'20.gmapp8 -O ./outputdump/$i.fasta

done
./create_index -G ./data/merged-file.fasta -I ./index/merged-file
./mappability -I ./index/merged-file -K $k -E $e -o $o -O ./output/merged-file
./mappability_dump -I ./output/merged-file_${e}_${k}_20.gmapp8 -O ./merged-file.fasta
