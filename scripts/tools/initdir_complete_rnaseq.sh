#!/usr/bin/env bash
# Ensures correct file hierarchy is present in data folder.
#INPUT: data folder directory

#ls interim/rna/
#aligned sorted_index mark_duplicates splitNtrim 

#ls processed/rna/
#recalibrated

set -u
set -o errexit
#set -o pipefail
set -o nounset
#set -o xtrace
#set -o verbose

DATA_DIRECTORY=/data/local/proj/bioinformatics_project/data

if [ ! -d "$DATA_DIRECTORY" ]
then
	echo "ERROR: Data directory not found"
	exit 1    
fi

echo "Data directory found"

RNA_INTERIM_DIR="$DATA_DIRECTORY"/interim/rna
if [ ! -d "$RNA_INTERIM_DIR/aligned" ]
then
	mkdir -p  "$RNA_INTERIM_DIR"/aligned

fi

if [ ! -d "$RNA_INTERIM_DIR"/sorted_index ]
then
	mkdir -p "$RNA_INTERIM_DIR/sorted_index"
fi

if [ ! -d "$RNA_INTERIM_DIR"/mark_duplicates ]
then
	mkdir -p "$RNA_INTERIM_DIR/mark_duplicates"
fi

if [ ! -d "$RNA_INTERIM_DIR"/splitNtrim ]
then
	mkdir -p "$RNA_INTERIM_DIR/splitNtrim"
fi


RNA_PROCESSED_DIR="$DATA_DIRECTORY"/processed/rna

if [ ! -d "$RNA_PROCESSED_DIR"/recalibration ]
then
	mkdir -p "$RNA_PROCESSED_DIR"/recalibrated
fi
