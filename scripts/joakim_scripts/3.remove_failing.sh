#!/usr/bin/env bash
if [ ! -d "logs" ]
then
    mkdir logs
fi
exec 1>logs/$(basename "$0").$(date +"%r%d%h%y" | awk '{print $1"_"$2}')."$RANDOM".log 2>&1

echo "$CONDA_PREFIX"

set -u
set -o errexit
set -o pipefail
set -o nounset
set -o xtrace
set -o verbose

vcf_gz="$1"
outname=$(dirname "$vcf_gz")/$(basename "$vcf_gz" ".vcf.gz").pass.vcf
cat <(bgzip -d -c "$vcf_gz" | grep "#") <(bgzip -d -c "$vcf_gz" | grep -v "#" | grep -w PASS) > "$outname"
#bgzip -f "$outname"
#tabix -f "$outname".gz
