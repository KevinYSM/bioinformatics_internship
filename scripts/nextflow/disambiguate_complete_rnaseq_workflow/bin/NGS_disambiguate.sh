#!/usr/bin/env bash


set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam_human="$1"
bam_mouse="$2"


sname=$(basename "$bam_human" | awk -F"." '{print $1}')

ngs_disambiguate -s "$sname" -a star "$bam_human" "$bam_mouse"
