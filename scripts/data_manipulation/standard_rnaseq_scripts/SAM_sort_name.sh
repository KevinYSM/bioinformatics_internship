#!/usr/bin/env bash

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose


output_dir=$2
samtools sort -n -m 10G -O BAM -o $2.sorted_by_name.bam "$1"
