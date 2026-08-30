#!/bin/bash
# Extract candidate MVA gene regions from the full VCF into a small subset VCF for VEP.
# PREPARING THE VCF subset file for https://jun2026.archive.ensembl.org/Homo_sapiens/Tools/VEP
# Select: 
#   Frequency data for co-located variants > gnomAD (exomes) allele frequencies
#   Frequency data for co-located variants > gnomAD (genomes) allele frequencies
#   Phenotype data and citations > Phenotypes

# output: https://jun2026.archive.ensembl.org/Homo_sapiens/Tools/VEP/Results?tl=byFzcDajq2wkhzw9-12021070
set -euo pipefail

VCF=/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/track1/targeted
BEDFILE="$OUTDIR/step1_genes_bedfile.bed"
OUTVCF="$OUTDIR/step2_MVA_subset.vcf"
OUTVCFGZ="$OUTVCF.gz"

# Make sure the index is fresh (fixes the earlier .tbi warning)
tabix -p vcf -f "$VCF"

echo "=== Candidate gene regions ==="
cat "$BEDFILE"

# Extract just these regions from the full VCF
bcftools view -R "$BEDFILE" -O z -o "$OUTVCFGZ" "$VCF"
tabix -p vcf -f "$OUTVCFGZ"

echo ""
echo "=== Subset created: $OUTVCFGZ ==="
echo "Variant count:"
bcftools view -H "$OUTVCFGZ" | wc -l

zcat $OUTVCFGZ > $OUTVCF