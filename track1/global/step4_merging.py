import pandas as pd

STEP1_OUT_GENOTYPE_FILE = "./step1_out.tsv"
STEP3_RANKED_CSV = "./step3_ranked_candidates.csv"
STEP4_MERGED_CSV = "./step4_merged.csv"

genotypes = pd.read_csv(STEP1_OUT_GENOTYPE_FILE, sep="\t", dtype={"CHROM": str}, low_memory=False)
candidates = pd.read_csv(STEP3_RANKED_CSV)

# Extract chromosome and position from Location
candidates[["CHROM", "POS"]] = candidates["Location"].str.extract(
    r"^([^:]+):(\d+)"
)
candidates["POS"] = candidates["POS"].astype(int)
genotypes["CHROM"] = genotypes["CHROM"].astype(str)

merged = candidates.merge(
    genotypes.rename(columns={"REF": "REF_ALLELE",}),
    on=["CHROM", "POS", "REF_ALLELE"],
    how="left"
)


# --- Simple hard filters on the new genotype columns ---
merged = merged[merged["FILTER"] == "PASS"] # drop calls that failed GATK's own filters
merged = merged[merged["GT"].isin(["0/1", "1/1", "1/0", "1|0", "0|1", "1|1"])] # keep only het/hom-alt calls (child actually carries the allele)
merged = merged[merged["ALT"].notna()] # drop rows with no alternate allele recorded
merged = merged[pd.to_numeric(merged["AD"].str.split(",").str[1], errors="coerce").fillna(0) > 0] # require at least 1 read supporting the alt allele
merged = merged[pd.to_numeric(merged["DP"], errors="coerce") >= 10]  # minimum total read depth to trust the call

merged.to_csv(STEP4_MERGED_CSV, index=False)
