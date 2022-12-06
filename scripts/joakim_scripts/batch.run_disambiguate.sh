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

while read line
do
    bam_human=$(echo "$line" | awk '{print $1}')
    bam_mouse=$(echo "$line" | awk '{print $2}')
    ./3.run_disambiguate.sh "$bam_human" "$bam_mouse" &
done < pairs.txt
wait
