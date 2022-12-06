#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}').log 2>&1

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

bam_human="$1"
bam_mouse="$2"

outdir_base=results/ngs_disambiguate
sname=$(basename "$bam_human" | awk -F"." '{print $1}')
outdir="$outdir_base"/"$sname"
if [ ! -d "$outdir" ]
then
    mkdir -p "$outdir"
fi
ngs_disambiguate -s "$sname" -o "$outdir" -a star "$bam_human" "$bam_mouse"
