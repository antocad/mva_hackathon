VCF=/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/track1/targeted
BEDFILE="$OUTDIR/step1_genes_bedfile.bed"
OUTFILE="$OUTDIR/step1_out.tsv"
# Candidate genes, GRCh38 no-chr coords (double check these against Ensembl/UCSC — I'm giving approximate ranges, verify exact exon boundaries)
# Get the spans here: https://genome.ucsc.edu/cgi-bin/hgSearch,  Step 1: Select a Genome: Human (hg38), Step 2: Search for a Term (see below)  (MANE Select)
# https://www.ncbi.nlm.nih.gov/gtr/tests/567860.2
# BUB1B:   15:40161069-40221123
# CEP57:   11:95790498-95832693
# TRIP13:   5:892884-918120
# https://www.ncbi.nlm.nih.gov/gtr/conditions/C4551972
# CENATAC: 11:118998138-119015793
# SLF2:    10:100912963-100965134
# SMC5:     9:70258978-70354873
# MAD1L1:   7:1815795-2232945

# BED file: chrom, start (0-based), end, gene name
cat > "$BEDFILE" << 'EOF'
15	40161068	40221123	BUB1B
11	95790497	95832693	CEP57 
5	892883	918120	TRIP13
11	118998137	119015793	CENATAC
10	100912962	100965134	SLF2
9	70258977	70354873	SMC5
7	1815794	2232945	MAD1L1
EOF

tabix -p vcf -f "$VCF"
# PASS vs filtered counts — gives you a sense of how aggressive the filtering was
# bcftools view -H "$VCF" | cut -f7 | sort | uniq -c | sort -rn

rm -f $OUTFILE
echo "CHROM\tPOS\tREF\tALT\tFILTER\tGT\tAD\tDP" > $OUTFILE
bcftools view -R "$BEDFILE" "$VCF" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\n' >> $OUTFILE