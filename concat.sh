
#! /usr/bin/bash
cd outputdump
for i in {0..31} ; do
sed 's/$/     /' $i.fasta > $i+1.fasta 
 
done 


cd outputdump
for i in {0..31} ; do
cat $i+1.fasta > all 
 
done 


exit 0


