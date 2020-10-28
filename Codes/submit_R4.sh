#!/bin/bash
files=`ls fileid.*`
for ii in $files;
do
echo $ii | cut -f2 -d"."
qsub submit_preR3.sh
sleep 2
done


