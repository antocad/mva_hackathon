VCF=/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz
OUTDIR=/home/molcure/storageB/antoine/mva_hackathon/track1/global
BEDFILE="$OUTDIR/step1_genes_bedfile.bed"
OUTFILE="$OUTDIR/step1_out.tsv"

tabix -p vcf -f "$VCF"
# PASS vs filtered counts — gives you a sense of how aggressive the filtering was
# bcftools view -H "$VCF" | cut -f7 | sort | uniq -c | sort -rn

rm -f $OUTFILE
echo "CHROM\tPOS\tREF\tALT\tFILTER\tGT\tAD\tDP" > $OUTFILE
bcftools view "$VCF" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\n' >> $OUTFILE