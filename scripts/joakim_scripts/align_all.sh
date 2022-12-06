#!/bin/bash

for f in $(find /data/local/data/GWA/something/*.fastq)
do
	align.sh $f & 
done
wait

