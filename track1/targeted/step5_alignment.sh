#!/bin/bash
set -euo pipefail

# Activate the conda environment once, for the whole script
source /home/molcure/miniconda3/etc/profile.d/conda.sh
conda activate antoine

FASTQ_DIR=/home/molcure/storageB/antoine/mva_hackathon/data/fastq
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/data/step5
GENOME_DIR="$OUTDIR/refGenome"

# # Download the genome 
# wget http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# mv Homo_sapiens.GRCh38.dna.primary_assembly.fa "$GENOME_DIR/step5_genome.fasta"

# # Index reference once
# bwa index "$GENOME_DIR/step5_genome.fasta"
# samtools faidx -@ 20 "$GENOME_DIR/step5_genome.fasta"
# picard CreateSequenceDictionary R="$GENOME_DIR/step5_genome.fasta" O="$GENOME_DIR/step5_genome.dict"

# # Align all 4 lanes, merging as you go (adjust sample name/read group as needed)
# for lane in L001 L002 L003 L004; do
#   bwa mem -t 20 -R "@RG\tID:${lane}\tSM:EX2312012\tPL:ILLUMINA" \
#     "${GENOME_DIR}/step5_genome.fasta" \
#     "${FASTQ_DIR}/WGS_EX2312012_HGWCNDSX7_S16_${lane}_R1_001.fastq.gz" \
#     "${FASTQ_DIR}/WGS_EX2312012_HGWCNDSX7_S16_${lane}_R2_001.fastq.gz" \
#     | samtools sort -@ 20 -o "${OUTDIR}/${lane}.sorted.bam" -
# done

# samtools merge -@ 20 $OUTDIR/merged.bam $OUTDIR/L001.sorted.bam $OUTDIR/L002.sorted.bam $OUTDIR/L003.sorted.bam $OUTDIR/L004.sorted.bam
# picard MarkDuplicates I=$OUTDIR/merged.bam O=$OUTDIR/final.bam M=$OUTDIR/dup_metrics.txt
# samtools index -@ 20 $OUTDIR/final.bam

samtools coverage $OUTDIR/final.bam > $OUTDIR/coverage_per_chrom.txt
cat $OUTDIR/coverage_per_chrom.txt

# cnvkit.py batch $OUTDIR/final.bam -m wgs \
#   -f $GENOME_DIR/step5_genome.fasta \
#   -n -p 20 \
#   --target-avg-size 10000 \
#   --output-reference $OUTDIR/flat_reference.cnn \
#   -d $OUTDIR/cnvkit_output/