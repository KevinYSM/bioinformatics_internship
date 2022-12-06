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

bam="$1"
#gtf=genome/Homo_sapiens.GRCh37.75.sorted.gtf
gtf=/data/local/reference/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf
htseq-count -r name -q -f bam -s yes -m intersection-strict "$bam" $gtf > "$bam".s_yes.gene_counts &
htseq-count -r name -q -f bam -s no -m intersection-strict "$bam" $gtf > "$bam".s_no.gene_counts &
htseq-count -r name -q -f bam -s reverse -m intersection-strict "$bam" $gtf > "$bam".s_reverse.gene_counts &
wait
