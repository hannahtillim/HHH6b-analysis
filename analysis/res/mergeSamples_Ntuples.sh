#!/bin/bash
# This script merges all the samples in a specified folder
# into combined background samples and NTuples.
rm -rf $1/background
mkdir $1/background

for d in $1*/
do

        for f in $d*NTuple.dat; do
			filename=$(basename $f)
        	if [ -f $1$filename ]; 
			then
            	echo "Total Ntuple merging " $f
            	cat $f | awk ' NR>1 {print;}' >> $1$filename 
            else
            	echo "Total NTuple generation " $f
            	cp $f $1$filename
        	fi
        done
done

