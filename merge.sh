#!/bin/bash

for i in {0..31}; do

cat ./data/"$i.fasta" >> ./data/merged-file.fasta

 done
