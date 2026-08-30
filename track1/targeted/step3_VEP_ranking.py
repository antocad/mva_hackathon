import sys
import re
import pandas as pd
import numpy as np

STEP2_VEP_OUTPUT_FILE = "/home/molcure/storageB/antoine/mva_hackathon/track1/targeted/step2_VEP_output.txt"

# Higher = more severe. Used to score whichever consequence(s) a row has.
CONSEQUENCE_RANK = {
    "stop_gained": 5,
    "frameshift_variant": 5,
    "splice_acceptor_variant": 4,
    "splice_donor_variant": 4,
    "start_lost": 4,
    "stop_lost": 4,
    "splice_region_variant": 2,
    "missense_variant": 1,
}
 
GNOMAD_AF_COLUMNS = ["gnomADe_AF", "gnomADg_AF", "AF"]  # AF = 1000 Genomes fallback
RARE_AF_THRESHOLD = 0.01  # above this = common, deprioritize
 
# Confirmed MVA candidate genes (add PLK4/CEP63/CEP192 as a tier-2 set if needed)
CANDIDATE_GENES = {"BUB1B", "CEP57", "TRIP13", "CENATAC", "SLF2", "SMC5", "MAD1L1"}

def normalize_text(s: str) -> str:
    return re.sub(r"[_\-]+", " ", str(s).lower()).strip()

def normalize_phenotype_text(series):
    return series.fillna("").str.lower().str.replace("_", " ", regex=False).str.replace("-", " ", regex=False)

PHENOTYPE_SYNONYMS = {
    "rhabdomyosarcoma":      ["rhabdomyosarcoma", "soft tissue sarcoma"],
    "nephrocalcinosis":      ["nephrocalcinosis", "renal calcification", "kidney calcification"],
    "short_stature":         ["short stature", "dwarfism", "growth retardation", "growth failure", "primordial dwarfism"],
    "failure_to_thrive":     ["failure to thrive", "poor weight gain", "growth failure", "faltering growth"],
    "muscle_atrophy":        ["muscle atrophy", "skeletal muscle atrophy", "muscle hypotonia", "low muscle mass"],
    "prematurity":           ["premature birth", "prematurity", "preterm birth"],
    "iugr":                  ["intrauterine growth restriction", "iugr", "fetal growth restriction",
                               "small for gestational age", "sga", "low birth weight"],
    "recurrent_miscarriage": ["recurrent miscarriage", "recurrent pregnancy loss",
                              "recurrent spontaneous abortion", "habitual abortion", "pregnancy loss"],
    "instability_disorders": ["aneuploidy", "mosaic variegated aneuploidy", "mva", "chromosome instability",
                               "chromosomal instability", "fanconi anemia", "bloom syndrome",
                               "nijmegen breakage", "cancer predisposition", "bone marrow failure",
                               "seckel", "cohesinopathy"],
}
PHENOTYPE_SYNONYMS = {k: [normalize_text(t) for t in v] for k, v in PHENOTYPE_SYNONYMS.items()}
PHENOTYPE_PATTERNS = {k: "|".join(v) for k, v in PHENOTYPE_SYNONYMS.items()}

PHENOTYPE_WEIGHTS = {
    "rhabdomyosarcoma":      2.0,   # the oncologic trigger feature itself
    "nephrocalcinosis":      1.0,   # renal anomaly
    "short_stature":         1.0,   # growth restriction
    "failure_to_thrive":     1.0,   # growth restriction
    "muscle_atrophy":        1.0,   # growth restriction
    "prematurity":           1.5,   # adverse perinatal history
    "iugr":                  1.5,   # adverse perinatal history
    "recurrent_miscarriage": 2.5,   # explicitly flagged as clinically significant, not incidental
    "instability_disorders": 1.5,   # mechanism-level corroboration
}

CLINVAR_POSITIVE_TERMS = {"pathogenic", "likely_pathogenic"}
CLINVAR_CONFLICTING_TERMS = {"conflicting_interpretations_of_pathogenicity", "conflicting_classifications_of_pathogenicity"}

def classify_clinvar(clin_sig):
    if not isinstance(clin_sig, str) or not clin_sig or clin_sig == "-":
        return "none"
    # VEP separates multiple ClinVar values with & or , or /
    terms = {t.strip().lower() for t in re.split(r"[,&/]", str(clin_sig))}
    if terms & CLINVAR_POSITIVE_TERMS:
        return "pathogenic"
    if terms & CLINVAR_CONFLICTING_TERMS:
        return "conflicting"
    return "other"

def has_clinvar(clin_sig):
    return classify_clinvar(clin_sig) == "pathogenic"

########################################################################
########################################################################
########################################################################
########################################################################

# Load VEP output txt file
with open(STEP2_VEP_OUTPUT_FILE, "r") as f:
    header = f.readline().lstrip("#").rstrip("\n").split("\t")

df = pd.read_csv(
    STEP2_VEP_OUTPUT_FILE,
    sep="\t",
    skiprows=1,
    names=header,
    low_memory=False
)

gene_col = "SYMBOL"
text_cols = [c for c in ["Consequence", "CLIN_SIG", "PHENOTYPES", "SIFT", "PolyPhen"] if c in df.columns]
df[text_cols] = df[text_cols].astype("string")

# --- gnomAD AF: max across whichever frequency columns exist, row-wise ---
af_cols = [c for c in GNOMAD_AF_COLUMNS if c in df.columns]
df["gnomAD_AF_max"] = df[af_cols].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)

# --- Consequence: VEP packs multiple values per row as "a,b,c" ---
# A single vectorized regex test tells us if ANY damaging term is present.
damaging_pattern = "|".join(CONSEQUENCE_RANK.keys())
df["is_damaging_consequence"] = df["Consequence"].str.contains(
    damaging_pattern, case=False, na=False
)

# Consequence severity score: take the max rank among any matching terms.
# Building one column per term and taking a row-wise max is fully vectorized.
rank_cols = pd.DataFrame(
    {
        term: df["Consequence"].str.contains(term, case=False, na=False) * rank
        for term, rank in CONSEQUENCE_RANK.items()
    }
)
df["consequence_score"] = rank_cols.max(axis=1)

# phenotype score
if "PHENOTYPES" in df.columns:
    df["PHENOTYPES_norm"] = normalize_phenotype_text(df["PHENOTYPES"])
    phenotype_cat_cols = pd.DataFrame({
        key: df["PHENOTYPES_norm"].str.contains(pattern, case=False, na=False, regex=True) * PHENOTYPE_WEIGHTS[key]
        for key, pattern in PHENOTYPE_PATTERNS.items()
    })
    df["phenotype_score"] = phenotype_cat_cols.sum(axis=1)
else:
    df["phenotype_score"] = 0.0

df["is_candidate_gene"] = df[gene_col].isin(CANDIDATE_GENES)

df["is_rare"] = df["gnomAD_AF_max"] < RARE_AF_THRESHOLD

clin_col = "CLIN_SIG" if "CLIN_SIG" in df.columns else None
df["has_clinvar_hit"] = (
    df[clin_col].apply(has_clinvar) if clin_col else False
)

# SIFT: "deleterious(0.01)" / PolyPhen: "probably_damaging(0.95)"
sift_bad = (
    df["SIFT"].str.contains("deleterious", case=False, na=False)
    if "SIFT" in df.columns
    else False
)
polyphen_bad = (
    df["PolyPhen"].str.contains("damaging", case=False, na=False)
    if "PolyPhen" in df.columns
    else False
)
df["predicted_damaging"] = sift_bad | polyphen_bad

# --- Filter to plausible candidates ---
filtered = df[df["is_damaging_consequence"] & df["is_rare"] & df["is_candidate_gene"]].copy()

# --- Priority score: ClinVar hit > SIFT/PolyPhen > known phenotype > consequence severity > rarity ---
filtered["priority_score"] = (
    filtered["has_clinvar_hit"].astype(int) * 100
    + filtered["predicted_damaging"].astype(int) * 15
    + filtered["phenotype_score"] * 10
    + filtered["consequence_score"] * 10
    + (1 - filtered["gnomAD_AF_max"]) * 5
)

ranked = filtered.sort_values("priority_score", ascending=False)

display_cols = [c for c in ["Uploaded_variation",
                            "Location",
                            gene_col,
                            "REF_ALLELE", # new
                            "UPLOADED_ALLELE", # new
                            "Consequence",
                            "gnomAD_AF_max",
                            "SIFT",
                            "PolyPhen",
                            clin_col,
                            "PHENOTYPES",
                            "priority_score"] if c and c in ranked.columns
]

final_ranked = ranked[display_cols]

sub_dfs = []
for loc in set(final_ranked.Location):
    sub_df = final_ranked[final_ranked.Location == loc]
    sub_df.Consequence = (",".join(set((",".join(sub_df.Consequence.to_numpy())).split(","))))
    sub_df.CLIN_SIG = (",".join(set((",".join(sub_df.CLIN_SIG.to_numpy())).split(","))))
    sub_df.PHENOTYPES=(",".join(set((",".join(sub_df.PHENOTYPES.to_numpy())).split(","))))
    sub_dfs.append(sub_df)

final_ranked = pd.concat(sub_dfs)
final_ranked.sort_values("priority_score", inplace=True, ascending=False)
final_ranked.drop_duplicates(subset=["Location", "SYMBOL"], inplace=True, keep="first")
final_ranked.reset_index(drop=True, inplace=True)

print(f"\n{len(final_ranked)} candidate variants after filtering (rare + damaging + in candidate gene list)")
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)
# print(final_ranked.to_string(index=False))

final_ranked.to_csv("step3_ranked_candidates.csv", index=False)
print("\nSaved full ranked list to ranked_candidates.csv")

