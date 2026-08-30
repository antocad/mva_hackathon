#!/bin/bash

# tmux new -s vep 'bash /home/molcure/storageB/antoine/mva_hackathon/track1/global/step2_VEP.sh; bash'

# Extract candidate MVA gene regions from the full VCF into a small subset VCF for VEP.
# PREPARING THE VCF subset file for https://jun2026.archive.ensembl.org/Homo_sapiens/Tools/VEP

set -euo pipefail

VCF=/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz
CACHE=/home/molcure/storageB/antoine/mva_hackathon/vep
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/data/step2_global/out
mkdir -p "$OUTDIR"

GTF_GZ=/home/molcure/storageB/antoine/mva_hackathon/data/step2_global/Homo_sapiens.GRCh38.110.gtf.gz
GENE_BED=/home/molcure/storageB/antoine/mva_hackathon/data/step2_global/genes_genome_wide_pad20.bed
OUTVCF="$OUTDIR/step2_global_subset.vcf.gz"

# Make sure the index is fresh
tabix -p vcf -f "$VCF"

echo "Total variant count (unfiltered):"
bcftools view -H "$VCF" | wc -l

# --- 1. Build a genome-wide "genic" BED: every gene body, ±20bp for splice sites ---
# Same Ensembl release/build as your MANE/VEP annotation, so coordinates stay consistent
# with the rest of your pipeline (GRCh38, no 'chr' prefix).
if [ ! -f "$GENE_BED" ]; then
  wget -q -O "$GTF_GZ" \
    "http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

  zcat "$GTF_GZ" \
    | awk 'BEGIN{OFS="\t"} $0 !~ /^#/ && $3=="gene" {
        start = $4 - 21; if (start < 0) start = 0;
        end = $5 + 20;
        print $1, start, end
      }' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -i - > "$GENE_BED"
fi

echo "Merged genic intervals:"
wc -l "$GENE_BED"

# --- 2. PASS-only + genic-region subset, in one bcftools pass ---
# (uncomment and run the lines below if you need to regenerate the subset VCF)
bcftools view -f PASS -R "$GENE_BED" "$VCF" -Oz -o "$OUTVCF"
tabix -p vcf -f "$OUTVCF"

# Pre-run checks: ensure the subset VCF exists and is indexed for VEP
if [ ! -f "$OUTVCF" ]; then
  echo "ERROR: required VCF not found: $OUTVCF" >&2
  echo "Please create the subset VCF before running VEP (see commented bcftools lines)." >&2
  exit 1
fi

# Ensure index (.tbi) exists; create if missing
if [ ! -f "${OUTVCF}.tbi" ]; then
  if command -v tabix >/dev/null 2>&1; then
    echo "Index missing for $OUTVCF — creating ${OUTVCF}.tbi"
    tabix -p vcf "$OUTVCF"
  else
    echo "ERROR: tabix not available to create index ${OUTVCF}.tbi" >&2
    echo "Please install samtools/tabix or create the index manually." >&2
    exit 1
  fi
fi

# 1. Get the VEP_plugins repo (one-time)
git clone https://github.com/Ensembl/VEP_plugins.git "$CACHE/Plugins"

# 2. Pre-download the phenotype annotation file — required because --offline
#    mode will refuse to auto-download it on first use (it dies rather than fetch)
mkdir -p "$CACHE/Plugins/phenotype_data"
docker run --user $(id -u):$(id -g) -v "$CACHE":/opt/vep/.vep \
  ensemblorg/ensembl-vep:release_110.1 \
  vep --species homo_sapiens --plugin Phenotypes,dir="/opt/vep/.vep/Plugins/phenotype_data" \
      --cache --dir_cache /opt/vep/.vep \
      -i /opt/vep/.vep/Plugins/dummy.vcf -o STDOUT

docker run --user $(id -u):$(id -g) -v "$CACHE":/opt/vep/.vep \
    -v "$OUTVCF":/indata/step2_global_subset.vcf.gz:ro \
    -v "${OUTVCF}.tbi":/indata/step2_global_subset.vcf.gz.tbi:ro \
    -v "$OUTDIR":/out \
  ensemblorg/ensembl-vep:release_110.1 \
  vep --input_file /indata/step2_global_subset.vcf.gz \
      --output_file /out/step2_VEP_output.txt \
      --tab \
      --offline --cache --dir_cache /opt/vep/.vep \
      --assembly GRCh38 --fork 16 --pick \
      --symbol --show_ref_allele --uploaded_allele \
      --af_gnomade --af_gnomadg --af \
      --sift b --polyphen b \
      --check_existing --clin_sig_allele 1 \
      --dir_plugins /opt/vep/.vep/Plugins \
      --plugin Phenotypes,dir=/opt/vep/.vep/Plugins/phenotype_data \
      --gene_phenotype \
      --fields "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,Gene,REF_ALLELE,UPLOADED_ALLELE,SIFT,PolyPhen,AF,gnomADe_AF,gnomADg_AF,CLIN_SIG,PHENOTYPES,Existing_variation"