"""
step5_build_submission.py

Builds the Track 1 submission CSV from step4_merged.csv (the ranked,
genotype-filtered candidate table produced by step3/step4).

Output columns (per hackathon spec):
    proband_id, chrom_1, pos_1, ref_1, alt_1,
    chrom_2, pos_2, ref_2, alt_2,
    epcr, finding_type, notes

- One row per proposed causal variant, OR one row per compound-het pair
  (chrom_2/pos_2/ref_2/alt_2 filled in only for pairs).
- Up to 10 rows total, ranked by EPCR (Estimated Probability of Causal
  Relationship), highest = most likely causal. EPCR only needs to be
  correctly *ordered* relative to your own other guesses -- there's no
  fixed "confident" threshold.
- finding_type: "primary" if the variant/pair sits in a known/curated
  MVA-relevant gene or already has a ClinVar pathogenic hit; otherwise
  "secondary" (incidental / exploratory finding from the genome-wide pass).
  Review and adjust these by hand before final submission -- this is a
  reasonable automatic default, not a clinical judgment.
"""

import re
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------
# Config -- EDIT THESE before running
# ----------------------------------------------------------------------
STEP4_FILES = {
    "targeted": "/home/molcure/storageB/antoine/mva_hackathon/track1/targeted/step4_merged.csv",
    "global":   "/home/molcure/storageB/antoine/mva_hackathon/track1/global/step4_merged.csv",
}

PROBAND_ID = "PROBAND01"   # from the hackathon dataset docs `bcftools query -l /home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz`

USERNAME = "antocad3"
APPROACH_NAME = "ensemble_without_alignment"
OUTPUT_FILE = f"{USERNAME}_{APPROACH_NAME}.csv"

MAX_ROWS = 10

# Known/curated MVA-relevant genes -- used only to set a sensible default
# for finding_type (primary vs secondary), not to filter candidates.
KNOWN_MVA_GENES = {
    "BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1",
}
# KNOWN_MVA_GENES = {
#     "BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1",
#     # Tier 2 differential-diagnosis panel
#     "FANCA", "FANCC", "FANCG", "FANCD2", "BRCA2", "PALB2", "BRIP1",
#     "FANCE", "FANCF", "FANCI", "FANCL", "FANCB",
#     "NBN", "MRE11", "RAD50", "LIG4", "ATM", "BLM",
#     # Tier 3 differential-diagnosis panel
#     "ESCO2", "NIPBL", "SMC1A", "SMC3", "RAD21", "HDAC8", "RECQL4", "WRN", "ATR", "ATRI", "PCNT", "CEP152"
# }

OVERLAP_BONUS_PER_EXTRA_METHOD = 50

# Phase-confirmation discount: since read-based/parental phasing is being
# skipped for time, every compound-het pair here is *unconfirmed* -- we know
# both variants exist and are rare+damaging, but not whether they sit in
# trans (a real compound-het) or cis (same copy, other copy still intact,
# not actually causal). Discount reflects that added uncertainty relative to
# a homozygous call, which needs no such assumption. Tune as needed.
PHASE_UNCONFIRMED_PENALTY = 0.85


# ----------------------------------------------------------------------
# Load + tag each source
# ----------------------------------------------------------------------
def parse_location(loc: str):
    m = re.match(r"^(\w+):(\d+)(?:-(\d+))?$", str(loc))
    if not m:
        return None, None
    chrom, start, _ = m.groups()
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return chrom, int(start)

def classify_zygosity(gt):
    if pd.isna(gt) or str(gt) in (".", "./.", ".|."):
        return "no_call"
    alleles = re.split(r"[/|]", str(gt))
    if len(set(alleles)) == 1:
        return "hom_ref" if alleles[0] == "0" else "hom_alt"
    return "het"


frames = []
for source_name, path in STEP4_FILES.items():
    d = pd.read_csv(path, low_memory=False)
    required = ["Location", "REF_ALLELE", "ALT", "SYMBOL", "priority_score", "GT"]
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise ValueError(f"{path} is missing expected column(s): {missing}")
    d[["chrom", "pos"]] = d["Location"].apply(lambda x: pd.Series(parse_location(x)))
    d = d.dropna(subset=["chrom", "pos"]).copy()
    d["pos"] = d["pos"].astype(int)
    d["ref"] = d["REF_ALLELE"]
    d["alt"] = d["ALT"]
    d["zygosity"] = d["GT"].apply(classify_zygosity)
    d = d[d["zygosity"].isin(["het", "hom_alt"])].copy()  # drop no_call / hom_ref rows
    d["source"] = source_name
    print(f"Loaded {len(d)} candidate rows from '{source_name}' ({path})")
    frames.append(d)
 
all_variants = pd.concat(frames, ignore_index=True)
all_variants["variant_key"] = (
    all_variants["chrom"] + ":" + all_variants["pos"].astype(str) + ":"
    + all_variants["ref"].astype(str) + ">" + all_variants["alt"].astype(str)
)


# ----------------------------------------------------------------------
# Collapse to one row per unique variant, tracking which method(s) found it
# ----------------------------------------------------------------------
def collapse_variant_group(grp: pd.DataFrame) -> pd.Series:
    methods = sorted(grp["source"].unique())
    n_methods = len(methods)
    best = grp.loc[grp["priority_score"].idxmax()]  # inherit annotation from best-scoring hit

    zygosities = set(grp["zygosity"].unique())
    zygosity_conflict = len(zygosities) > 1

    consequence = (",".join(set((",".join(grp.Consequence.to_numpy())).split(","))))
    clin_sig = (",".join(set((",".join(grp.CLIN_SIG.to_numpy())).split(","))))
    phenotypes = (",".join(set((",".join(grp.PHENOTYPES.str.lower().to_numpy())).split(","))))

    return pd.Series({
        "chrom": best["chrom"], "pos": best["pos"], "ref": best["ref"], "alt": best["alt"],
        "gene": best["SYMBOL"],
        "zygosity": best["zygosity"],
        "zygosity_conflict": zygosity_conflict,
        "n_methods": n_methods,
        "methods": "+".join(methods),
        "consequence": consequence,
        "gnomad_af": best.get("gnomAD_AF_max", np.nan),
        "clin_sig": clin_sig,
        "phenotypes": phenotypes,
        "combined_score": grp["priority_score"].sum(), # The sum of priority scores (if it appears in both methods then it has strong score)
    })


variants = (
    all_variants.groupby("variant_key", group_keys=False)
    .apply(collapse_variant_group)
    .sort_values(by="combined_score", ascending=False)
    .reset_index(drop=True)
)
print(f"\n{len(variants)} unique variants after merging both methods "
      f"({(variants['n_methods'] > 1).sum()} found by both)")

def is_known_gene(symbol):
    return str(symbol).upper() in KNOWN_MVA_GENES


CLINVAR_POSITIVE_TERMS = {"pathogenic", "likely_pathogenic"}
CLINVAR_CONFLICTING_TERMS = {"conflicting_interpretations_of_pathogenicity", "conflicting_classifications_of_pathogenicity"}

def classify_clinvar(clin_sig):
    if not isinstance(clin_sig, str) or not clin_sig or clin_sig == "-":
        return "none"
    # VEP separates multiple ClinVar values with & or ,
    terms = {t.strip().lower() for t in re.split(r"[,&/]", clin_sig)}
    if terms & CLINVAR_POSITIVE_TERMS:
        return "pathogenic"
    if terms & CLINVAR_CONFLICTING_TERMS:
        return "conflicting"
    return "other"

def has_clinvar(clin_sig):
    return classify_clinvar(clin_sig) == "pathogenic"

def build_notes(row, extra=""):
    parts = [extra] if extra else []
    parts.append(f"found_by={row['methods']}")
    if row["zygosity_conflict"]:
        parts.append("WARNING: zygosity call disagrees between methods -- verify manually")
    if pd.notna(row.get("consequence")) and row["consequence"]:
        parts.append(f"consequence={row['consequence']}")
    if pd.notna(row.get("gnomad_af")):
        parts.append(f"gnomAD_AF={row['gnomad_af']:.2e}")
    if isinstance(row.get("clin_sig"), str) and row["clin_sig"] and row["clin_sig"] != "nan":
        parts.append(f"ClinVar={row['clin_sig']}")
    if pd.notna(row.get("phenotypes")):
        parts.append(f"phenotypes={row['phenotypes']}")
    return "; ".join(parts)

# ----------------------------------------------------------------------
# 1. Homozygous single-variant candidates
# ----------------------------------------------------------------------
hom = variants[variants["zygosity"] == "hom_alt"].copy()
hom_rows = []
for _, r in hom.iterrows():
    hom_rows.append({
        "score": r["combined_score"],
        "chrom_1": r["chrom"], "pos_1": r["pos"], "ref_1": r["ref"], "alt_1": r["alt"],
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "is_mva_gene": is_known_gene(r["gene"]),
        "has_clinvar": has_clinvar(r["clin_sig"]),
        "notes": "gene: " + r['gene'] + " | 1st variant: " +build_notes(r),
    })


# ----------------------------------------------------------------------
# 2. Compound-het candidate pairs (within each gene, top 2 by combined_score)
# ----------------------------------------------------------------------
het = variants[variants["zygosity"] == "het"].copy()
# Treat blanks/placeholders as missing, not as a real shared "gene"
het["gene"] = het["gene"].replace({"": np.nan, "-": np.nan, "None": np.nan})
het_known = het.dropna(subset=["gene"])
het_unknown = het[het["gene"].isna()]

pair_rows = []
paired_indices = set()
for gene, grp in het_known.groupby("gene"):
    if gene == "MAD1L1" or gene == "BUB1B":
        print(grp)
    if len(grp) < 2:
        continue
    top2 = grp.sort_values("combined_score", ascending=False).head(2)
    v1, v2 = top2.iloc[0], top2.iloc[1]
    paired_indices.update(top2.index)
    pair_rows.append({
        "score": v1["combined_score"] + v2["combined_score"],
        "chrom_1": v1["chrom"], "pos_1": v1["pos"], "ref_1": v1["ref"], "alt_1": v1["alt"],
        "chrom_2": v2["chrom"], "pos_2": v2["pos"], "ref_2": v2["ref"], "alt_2": v2["alt"],
        "is_mva_gene": is_known_gene(gene),
        "has_clinvar": has_clinvar(v1["clin_sig"]) or has_clinvar(v2["clin_sig"]),
        "notes": (
            "gene: " + gene + " | 1st variant: " + build_notes(v1) + " | 2nd variant: " + build_notes(v2)
        ),
    })

# Any heterozygous variant that wasn't part of a pair -- either its gene had
# no second hit, or the gene symbol itself was unknown -- still gets a
# standalone row rather than being dropped. It doesn't fit a clean
# homozygous/compound-het recessive model on its own, so that caveat is
# flagged explicitly in the notes for manual review.
het_unpaired = pd.concat([
    het_known[~het_known.index.isin(paired_indices)],
    het_unknown,
])
for _, r in het_unpaired.iterrows():
    reason = "no gene symbol assigned (intergenic/unannotated)" if pd.isna(r["gene"]) else \
             "single heterozygous hit -- no second damaging variant found in this gene"
    pair_rows.append({
        "score": r["combined_score"],
        "chrom_1": r["chrom"], "pos_1": r["pos"], "ref_1": r["ref"], "alt_1": r["alt"],
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "is_mva_gene": is_known_gene(r["gene"]) if pd.notna(r["gene"]) else False,
        "has_clinvar": has_clinvar(r["clin_sig"]),
        "notes": "gene: " + (r['gene'] if pd.notna(r['gene']) else "<unknown>") + " | 1st variant: " + build_notes(r, extra=f"single heterozygous ({reason}; "
                                       f"does not independently fit the family's "
                                       f"apparent recessive inheritance pattern)"),
    })


# ----------------------------------------------------------------------
# Combine, rank, assign EPCR + finding_type, truncate to MAX_ROWS
# ----------------------------------------------------------------------
candidates = pd.DataFrame(hom_rows + pair_rows)
if candidates.empty:
    raise ValueError("No homozygous or compound-het candidates survived filtering.")

candidates = candidates.sort_values("score", ascending=False).reset_index(drop=True)
candidates = candidates.head(MAX_ROWS).copy()

min_s, max_s = candidates["score"].min(), candidates["score"].max()
if max_s > min_s:
    candidates["epcr"] = 0.05 + 0.95 * (candidates["score"] - min_s) / (max_s - min_s)
else:
    candidates["epcr"] = 1.0

candidates["finding_type"] = np.where(
    candidates["is_mva_gene"], "primary", "secondary"
)
candidates["proband_id"] = PROBAND_ID

output_cols = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "epcr", "finding_type", "notes",
]
submission = candidates[output_cols]

submission.to_csv(OUTPUT_FILE, index=False)
print(f"\nWrote {len(submission)} ranked candidate row(s) to {OUTPUT_FILE}")
print(submission.to_string(index=False))
