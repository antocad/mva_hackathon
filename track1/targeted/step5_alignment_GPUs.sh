#!/bin/bash
set -euo pipefail
 
# =========================================================================
# Step 5 (rewritten): GPU-accelerated alignment + sort + markdup
# via NVIDIA Parabricks fq2bam, replacing the CPU bwa/samtools/picard chain.
# =========================================================================
 
FASTQ_DIR=/home/molcure/storageB/antoine/mva_hackathon/data/fastq
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/data/step5
GENOME_DIR="$OUTDIR/refGenome"
PARABRICKS_IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"   # check NGC catalog for the current latest tag
 
mkdir -p "$OUTDIR" "$GENOME_DIR"
 
LOG="$OUTDIR/step5_$(date +%Y%m%d_%H%M%S).log"
echo "Logging to $LOG"
exec > >(tee -a "$LOG") 2>&1


# -------------------------------------------------------------------------
# 0. Reference genome: reuse an existing indexed FASTA if you already have
#    one from earlier in the pipeline (e.g. from VEP/annotation steps).
#    Only download if nothing is found.
# -------------------------------------------------------------------------
EXISTING_FASTA=$(find /home/molcure/storageB -iname "*GRCh38*primary_assembly*.fa" 2>/dev/null | head -n 1 || true)
 
if [[ -n "$EXISTING_FASTA" ]]; then
    echo "Found existing reference FASTA: $EXISTING_FASTA"
    REF_FASTA="$EXISTING_FASTA"
else
    echo "No existing reference found — downloading Ensembl GRCh38 release 110"
    wget -q --show-progress -O "$GENOME_DIR/step5_genome.fa.gz" \
        http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    gunzip -f "$GENOME_DIR/step5_genome.fa.gz"
    REF_FASTA="$GENOME_DIR/step5_genome.fasta"
    mv "$GENOME_DIR/step5_genome.fa" "$REF_FASTA" 2>/dev/null || true
fi

# BWA index files are still required by fq2bam (it uses bwa-mem2 under the hood
# but expects standard index files alongside the FASTA). Build them if missing.
if [[ ! -f "${REF_FASTA}.bwt" && ! -f "${REF_FASTA}.bwt.2bit.64" ]]; then
    echo "Building BWA index (one-time cost, ~1-2h on CPU)..."
    source /home/molcure/miniconda3/etc/profile.d/conda.sh
    conda activate antoine
    bwa index "$REF_FASTA"
fi
if [[ ! -f "${REF_FASTA}.fai" ]]; then
    source /home/molcure/miniconda3/etc/profile.d/conda.sh
    conda activate antoine
    samtools faidx "$REF_FASTA"
fi

# -------------------------------------------------------------------------
# 1. Confirm GPUs + Docker are visible before committing to the run
# -------------------------------------------------------------------------
echo "Checking GPU visibility..."
nvidia-smi --query-gpu=index,name,memory.total --format=csv

echo "Checking Docker GPU passthrough..."
docker run --rm --gpus all nvidia/cuda:12.2.0-base-ubuntu22.04 nvidia-smi --query-gpu=index,name --format=csv

# -------------------------------------------------------------------------
# 2. Pull Parabricks image (no-op if already pulled)
# -------------------------------------------------------------------------
echo "Pulling Parabricks image: $PARABRICKS_IMAGE"
docker pull "$PARABRICKS_IMAGE"

# -------------------------------------------------------------------------
# 3. Run fq2bam: does alignment + sorting + duplicate marking + indexing
#    in a single GPU-accelerated pass across all 4 lanes.
# -------------------------------------------------------------------------
REF_DIR=$(dirname "$REF_FASTA")
REF_NAME=$(basename "$REF_FASTA")

echo "Running Parabricks fq2bam on 2 GPUs..."
docker run --rm --gpus all \
    -v "$FASTQ_DIR":/fastq \
    -v "$REF_DIR":/ref \
    -v "$OUTDIR":/out \
    "$PARABRICKS_IMAGE" \
    pbrun fq2bam \
        --ref "/ref/${REF_NAME}" \
        --in-fq /fastq/WGS_EX2312012_HGWCNDSX7_S16_L001_R1_001.fastq.gz /fastq/WGS_EX2312012_HGWCNDSX7_S16_L001_R2_001.fastq.gz "@RG\tID:L001\tSM:EX2312012\tPL:ILLUMINA" \
        --in-fq /fastq/WGS_EX2312012_HGWCNDSX7_S16_L002_R1_001.fastq.gz /fastq/WGS_EX2312012_HGWCNDSX7_S16_L002_R2_001.fastq.gz "@RG\tID:L002\tSM:EX2312012\tPL:ILLUMINA" \
        --in-fq /fastq/WGS_EX2312012_HGWCNDSX7_S16_L003_R1_001.fastq.gz /fastq/WGS_EX2312012_HGWCNDSX7_S16_L003_R2_001.fastq.gz "@RG\tID:L003\tSM:EX2312012\tPL:ILLUMINA" \
        --in-fq /fastq/WGS_EX2312012_HGWCNDSX7_S16_L004_R1_001.fastq.gz /fastq/WGS_EX2312012_HGWCNDSX7_S16_L004_R2_001.fastq.gz "@RG\tID:L004\tSM:EX2312012\tPL:ILLUMINA" \
        --out-bam /out/final.bam \
        --num-gpus 2

echo "Done. Output: $OUTDIR/final.bam"
echo "final.bam is already coordinate-sorted, duplicate-marked, and indexed"
echo "(no separate samtools merge / picard MarkDuplicates / samtools index needed)."

ls -lh "$OUTDIR/final.bam" "$OUTDIR/final.bam.bai" 2>/dev/null || true