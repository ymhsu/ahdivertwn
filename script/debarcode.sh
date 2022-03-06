#!/bin/bash
file1="../data/seq/GA408/seq_list"

#debarcode using stacks
for n in $(cat "$file1")
do

first=$(echo $n | cut -d ":" -f 1)
second=$(echo $n | cut -d ":" -f 2)
third=$(echo $n | cut -d ":" -f 3)

echo ${first}
echo ${second}
echo ${third}

sudo process_radtags -p ${first} -b ${second} -o ${third} -e pstI -r -c -q -s 20 -i gzfastq

done




#sudo process_radtags -p ./20151013_HSMSR015/Raw/Sample_PCR_1-2_2/ -b ./barcode_lane1 -o ./debarcoded_lane1/first/ -e pstI -r -c -q -s 20 -i gzfastq
