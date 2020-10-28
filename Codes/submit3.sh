#!/bin/bash
i=2
while [ $i -lt 1002 ]
do
echo $i
aaa=$(sed -n "$i"'p' allgene.matrix.txt)
genename=$(echo $aaa | cut -d " " -f2)
genono=$(echo $aaa | cut -d " " -f3)
methylno=$(echo $aaa | cut -d " " -f4)

j=2
while [ $j -lt $((genono+1+1)) ]
do
k=$((genono+1+1))
while [ $k -lt $((genono+methylno+1+1)) ]
do
sed -n '1p' geneid_$((i-1))_$genename.txt>fileid.$((i-1)).$genename.$((j-1)).$((k-genono-1)).txt
sed -n "$j"'p' geneid_$((i-1))_$genename.txt>>fileid.$((i-1)).$genename.$((j-1)).$((k-genono-1)).txt
sed -n "$k"'p' geneid_$((i-1))_$genename.txt>>fileid.$((i-1)).$genename.$((j-1)).$((k-genono-1)).txt

#echo fileid.$((i-1)).$genename.$((j-1)).$((k-genono-1)).txt
#qsub submit_preR3.sh
#sleep 1.5
k=`expr $k + 1`
#Rscript add.R $genename.$j.$k.txt
done
j=`expr $j + 1`
done
i=`expr $i + 1`
done
