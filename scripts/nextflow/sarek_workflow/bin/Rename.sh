#!/usr/bin/env bash
for i in {1..22} X Y MT
do
echo "chr$i $i" > chr_name_conv.txt
        bcftools annotate --rename-chrs chr_name_conv.txt Schrodi_IL23_IL17_combined_RECAL_SNP_INDEL_variants.VA.chr$i.vcf.gz -Oz -o Schrodi_IL23_IL17_combined_RECAL_SNP_INDEL_variants.VA.chr$i.Minimac4.vcf.gz
        tabix -p vcf Schrodi_IL23_IL17_combined_RECAL_SNP_INDEL_variants.VA.chr$i.Minimac4.vcf.gz
done