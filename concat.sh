#!/bin/bash

for i in {0..31}; do


cat ./outputdump/"$i.fasta"  >> concatenated.fasta


 done
