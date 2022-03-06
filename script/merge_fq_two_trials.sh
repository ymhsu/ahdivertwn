#!/bin/bash

file1="../data/seq/GA408/debarcoded_lane1/lane1_list"


for lane1 in $(cat "$file1")
do

echo ${lane1}
gunzip -c ../data/seq/GA408/debarcoded_lane1/first/${lane1}.fq.gz >> ../data/seq/debarcoded_final/${lane1}_merge.fq
gunzip -c ../data/seq/GA408/debarcoded_lane1/second/${lane1}.fq.gz >> ../data/seq/debarcoded_final/${lane1}_merge.fq

done

file2="../data/seq/GA408/debarcoded_lane2/lane2_list"

for lane2 in $(cat "$file2")
do

echo ${lane2}
gunzip -c ../data/seq/GA408/debarcoded_lane2/first/${lane2}.fq.gz >> ../data/seq/debarcoded_final/${lane2}_merge.fq
gunzip -c ../data/seq/GA408/debarcoded_lane2/second/${lane2}.fq.gz >> ../data/seq/debarcoded_final/${lane2}_merge.fq

done

echo ${gz}
gzip ../data/seq/debarcoded_final/*.fq
