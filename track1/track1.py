#!/usr/bin/env python3
"""
track1.py

Unified Track 1 Analysis Pipeline for MVA Hackathon.
Executes genotype extraction, dockerized Ensembl VEP annotation, candidate ranking,
and genotype merging for both targeted (prior candidate genes) and global (genome-wide) approaches.
"""

import os
import sys
import re
import gzip
import argparse
import subprocess
from typing import Union, Dict, Tuple, Optional, Any, List
import pandas as pd
import numpy as np

# ----------------------------------------------------------------------
# Default priors: 7 MVA candidate genes with GRCh38 coordinates
# ----------------------------------------------------------------------
DEFAULT_MVA_GENES = {
    "BUB1B": "15:40161068-40221123",
    "CEP57": "11:95790497-95832693",
    "TRIP13": "5:892883-918120",
    "CENATAC": "11:118998137-119015793",
    "SLF2": "10:100912962-100965134",
    "SMC5": "9:70258977-70354873",
    "MAD1L1": "7:1815794-2232945",
}

# Consequence severity ranking
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

PHENOTYPE_WEIGHTS = {
    "rhabdomyosarcoma":      2.0,
    "nephrocalcinosis":      1.0,
    "short_stature":         1.0,
    "failure_to_thrive":     1.0,
    "muscle_atrophy":        1.0,
    "prematurity":           1.5,
    "iugr":                  1.5,
    "recurrent_miscarriage": 2.5,
    "instability_disorders": 1.5,
}

GNOMAD_AF_COLUMNS = ["gnomADe_AF", "gnomADg_AF", "AF"]
RARE_AF_THRESHOLD = 0.01

CLINVAR_POSITIVE_TERMS = {"pathogenic", "likely_pathogenic"}
CLINVAR_CONFLICTING_TERMS = {"conflicting_interpretations_of_pathogenicity", "conflicting_classifications_of_pathogenicity"}

def setup_cache(cache_dir: str, vep_image: str = "ensemblorg/ensembl-vep:release_110.1") -> None:
    """
    Prepare `cache_dir` with required VEP resources:
      - install the offline cache (with SIFT/PolyPhen) via VEP's INSTALL.pl inside Docker
      - clone VEP plugins into cache_dir/Plugins if missing
      - pre-fetch the Phenotypes plugin's GVF data (required before --offline runs)
    This is a no-op for resources that already exist.
    """
    cache_dir = os.path.abspath(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    user_id = f"{os.getuid()}:{os.getgid()}"

    # Offline cache (with SIFT/PolyPhen) is ~20GB. Must exist locally before the
    # Phenotypes pre-fetch below, otherwise VEP falls back to querying the public
    # Ensembl DB over the network, which can hang indefinitely.
    species_cache_dir = os.path.join(cache_dir, "homo_sapiens", "110_GRCh38")
    cache_ready = os.path.isdir(species_cache_dir) and bool(os.listdir(species_cache_dir))
    if not cache_ready:
        try:
            print(f"[Setup] Installing VEP offline cache (~20GB, may take a long time)...")
            installer_cmd = [
                "docker", "run", "--rm", "--user", user_id,
                "-v", f"{cache_dir}:/opt/vep/.vep",
                vep_image,
                "bash", "-lc",
                (
                    # INSTALL.pl lives in the image's source dir, not the cache dir, and must be run from there.
                    "cd /opt/vep/src/ensembl-vep && "
                    "perl INSTALL.pl -a cf "
                    "-s homo_sapiens "
                    "-y GRCh38 "
                    "-r 110 "
                    "--CACHEDIR /opt/vep/.vep "
                    "--NO_TEST"
                )
            ]
            subprocess.run(installer_cmd, check=True)
            cache_ready = os.path.isdir(species_cache_dir) and bool(os.listdir(species_cache_dir))
            print(f"[Setup] VEP cache installed at {species_cache_dir}." if cache_ready else f"[Setup] Warning: installer finished but {species_cache_dir} is still empty.")
        except FileNotFoundError:
            print(f"[Setup] Warning: docker not available; please run the VEP installer manually inside the VEP container.")
        except subprocess.CalledProcessError:
            print(f"[Setup] Warning: VEP cache installer failed. You may run the interactive installer inside the container to populate the cache.")
    else:
        print(f"[Setup] VEP cache already present at {species_cache_dir}. Skipping installer.")

    # Ensure VEP_plugins are present (clone only if missing)
    plugins_dir = os.path.join(cache_dir, "Plugins")
    if not os.path.exists(plugins_dir):
        try:
            print(f"[Setup] Cloning VEP plugins into {plugins_dir}...")
            subprocess.run(["git", "clone", "https://github.com/Ensembl/VEP_plugins.git", plugins_dir], check=True)
        except FileNotFoundError:
            print(f"[Setup] Warning: git not found. Please install git and clone VEP_plugins into {plugins_dir}.")
        except Exception:
            print(f"[Setup] Warning: failed to clone VEP plugins. You may clone manually into {plugins_dir}.")
    else:
        print(f"[Setup] VEP plugins already present in {plugins_dir}. Skipping clone.")

    # Pre-fetch the Phenotypes plugin's GVF data: in --offline mode VEP refuses to
    # auto-download it on first use, so it must be fetched with a non-offline dummy run.
    # Requires cache_ready, otherwise this falls back to the (hanging) public Ensembl DB.
    if not cache_ready:
        print(f"[Setup] Skipping Phenotypes plugin pre-fetch: offline cache not ready.")
        return

    phenotype_dir = os.path.join(plugins_dir, "phenotype_data")
    os.makedirs(phenotype_dir, exist_ok=True)
    phenotype_gvf = os.path.join(phenotype_dir, "Phenotypes.pm_homo_sapiens_110_GRCh38.gvf.gz")
    if not os.path.exists(phenotype_gvf):
        dummy_vcf = os.path.join(plugins_dir, "dummy.vcf")
        if not os.path.exists(dummy_vcf):
            with open(dummy_vcf, "w") as f:
                f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t100\t.\tA\tT\t.\t.\t.\n")
        try:
            print(f"[Setup] Pre-fetching Phenotypes plugin data (requires network)...")
            subprocess.run([
                "docker", "run", "--rm", "--user", user_id,
                "-v", f"{cache_dir}:/opt/vep/.vep",
                vep_image,
                "vep", "--species", "homo_sapiens",
                "--plugin", "Phenotypes,dir=/opt/vep/.vep/Plugins/phenotype_data",
                "--cache", "--dir_cache", "/opt/vep/.vep",
                "-i", "/opt/vep/.vep/Plugins/dummy.vcf", "-o", "STDOUT",
            ], check=True, capture_output=True)
            print(f"[Setup] Phenotypes plugin data ready.")
        except FileNotFoundError:
            print(f"[Setup] Warning: docker not available; please pre-fetch the Phenotypes plugin data manually.")
        except subprocess.CalledProcessError:
            print(f"[Setup] Warning: Phenotypes plugin data pre-fetch failed. Manual run may be required (see step2_VEP.sh).")
    else:
        print(f"[Setup] Phenotypes plugin data already present in {phenotype_dir}. Skipping download.")

def classify_clinvar(clin_sig: Any) -> str:
    if not isinstance(clin_sig, str) or not clin_sig or clin_sig == "-":
        return "none"
    terms = {t.strip().lower() for t in re.split(r"[,&/]", str(clin_sig))}
    if terms & CLINVAR_POSITIVE_TERMS:
        return "pathogenic"
    if terms & CLINVAR_CONFLICTING_TERMS:
        return "conflicting"
    return "other"

def has_clinvar(clin_sig: Any) -> bool:
    return classify_clinvar(clin_sig) == "pathogenic"

def normalize_text(s: str) -> str:
    return re.sub(r"[_\-]+", " ", str(s).lower()).strip()

def normalize_phenotype_text(series: pd.Series) -> pd.Series:
    return series.fillna("").str.lower().str.replace("_", " ", regex=False).str.replace("-", " ", regex=False)

PHENOTYPE_SYNONYMS = {k: [normalize_text(t) for t in v] for k, v in PHENOTYPE_SYNONYMS.items()}
PHENOTYPE_PATTERNS = {k: "|".join(v) for k, v in PHENOTYPE_SYNONYMS.items()}

# ----------------------------------------------------------------------
# Helper: Parse gene specifications
# ----------------------------------------------------------------------
def parse_genes_input(genes: Any) -> Optional[Dict[str, Tuple[str, int, int]]]:
    """
    Parses genes parameter into a dictionary: {gene_name: (chrom, start, end)}.
    Accepts:
      - None
      - Dict: {"BUB1B": "15:40161068-40221123"} or {"BUB1B": ("15", 40161068, 40221123)}
      - String: "BUB1B:15:40161068-40221123, CEP57:11:95790497-95832693" or "BUB1B, CEP57"
      - List/Set: ["BUB1B", "CEP57"] or [("BUB1B", "15:40161068-40221123")]
    """
    if genes is None:
        return None
    if isinstance(genes, str) and genes.strip().lower() in ("none", "null", ""):
        return None
    if isinstance(genes, dict):
        res = {}
        for g, val in genes.items():
            g_name = str(g).strip().upper()
            if isinstance(val, (tuple, list)):
                res[g_name] = (str(val[0]), int(val[1]), int(val[2]))
            elif isinstance(val, str):
                m = re.match(r"^(\w+)[:\s]+(\d+)[\-\s]+(\d+)$", val.strip())
                if m:
                    res[g_name] = (m.group(1), int(m.group(2)), int(m.group(3)))
                else:
                    raise ValueError(f"Cannot parse genomic range '{val}' for gene '{g}'")
            else:
                raise TypeError(f"Invalid range specification type for gene '{g}': {type(val)}")
        return res
    if isinstance(genes, str):
        res = {}
        items = re.split(r"[,;\n]+", genes)
        for item in items:
            item = item.strip()
            if not item:
                continue
            m = re.match(r"^([A-Za-z0-9_\-]+)[:\s]+(\w+)[:\s]+(\d+)[\-\s]+(\d+)$", item)
            if m:
                res[m.group(1).upper()] = (m.group(2), int(m.group(3)), int(m.group(4)))
            else:
                # Check single gene name in DEFAULT_MVA_GENES
                g_upper = item.split()[0].upper()
                if g_upper in DEFAULT_MVA_GENES:
                    sub_res = parse_genes_input({g_upper: DEFAULT_MVA_GENES[g_upper]})
                    if sub_res:
                        res.update(sub_res)
                else:
                    raise ValueError(f"Could not parse gene range specification: '{item}'")
        return res
    if isinstance(genes, (list, tuple, set)):
        res = {}
        for item in genes:
            sub = parse_genes_input(item)
            if sub:
                res.update(sub)
        return res
    raise TypeError(f"Unsupported type for genes argument: {type(genes)}")

def write_bed_file(parsed_genes: Dict[str, Tuple[str, int, int]], bed_path: str) -> str:
    with open(bed_path, "w") as f:
        for gene, (chrom, start, end) in parsed_genes.items():
            f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
    return bed_path

# ----------------------------------------------------------------------
# Step 1: Extract Genotypes
# ----------------------------------------------------------------------
def get_genotypes(
    vcf_path: str,
    out_dir: str,
    genes: Any = None,
    skip_if_exists: bool = False
) -> str:
    """
    Extracts genotype information from VCF into step1_out.tsv.
    If genes is provided, uses -R with a generated BED file.
    """
    os.makedirs(out_dir, exist_ok=True)
    outfile = os.path.join(out_dir, "step1_out.tsv")
    if skip_if_exists and os.path.exists(outfile) and os.path.getsize(outfile) > 0:
        print(f"[Step 1] Genotype output already exists at {outfile}, skipping extraction.")
        return outfile

    parsed_genes = parse_genes_input(genes)

    # Index VCF only if index does not exist
    tbi_path = vcf_path + ".tbi"
    if not os.path.exists(tbi_path):
        print(f"[Step 1] Creating tabix index for {vcf_path}...")
        subprocess.run(["tabix", "-p", "vcf", vcf_path], check=False)

    if parsed_genes is not None:
        bed_path = os.path.join(out_dir, "step1_genes_bedfile.bed")
        write_bed_file(parsed_genes, bed_path)
        cmd = (
            f"echo 'CHROM\tPOS\tREF\tALT\tFILTER\tGT\tAD\tDP' > '{outfile}' && "
            f"bcftools view -R '{bed_path}' '{vcf_path}' | "
            f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\n' >> '{outfile}'"
        )
    else:
        cmd = (
            f"echo 'CHROM\tPOS\tREF\tALT\tFILTER\tGT\tAD\tDP' > '{outfile}' && "
            f"bcftools view '{vcf_path}' | "
            f"bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD\t%DP]\n' >> '{outfile}'"
        )

    print(f"[Step 1] Extracting genotypes to {outfile}...")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    print(f"[Step 1] Genotypes extracted successfully.")
    return outfile

# ----------------------------------------------------------------------
# Step 2: Run Ensembl VEP (Dockerized)
# ----------------------------------------------------------------------
def run_vep(
    vcf_path: str,
    cache_dir: str,
    out_dir: str,
    method: str = "targeted",
    genes: Any = None,
    skip_if_exists: bool = False
) -> str:
    """
    Prepares subset VCF and executes dockerized Ensembl VEP.
    Outputs step2_VEP_output.txt.
    """
    os.makedirs(out_dir, exist_ok=True)
    vep_output = os.path.join(out_dir, "step2_VEP_output.txt")
    if skip_if_exists and os.path.exists(vep_output) and os.path.getsize(vep_output) > 0:
        print(f"[Step 2] VEP output already exists at {vep_output}, skipping execution.")
        return vep_output

    os.makedirs(cache_dir, exist_ok=True)
    plugins_dir = os.path.join(cache_dir, "Plugins")
    if not os.path.exists(plugins_dir):
        print(f"[Step 2] Cloning VEP plugins repository into cache dir...")
        subprocess.run(["git", "clone", "https://github.com/Ensembl/VEP_plugins.git", plugins_dir], check=True)

    phenotype_dir = os.path.join(plugins_dir, "phenotype_data")
    os.makedirs(phenotype_dir, exist_ok=True)

    parsed_genes = parse_genes_input(genes)

    if method == "targeted" or parsed_genes is not None:
        if parsed_genes is None:
            parsed_genes = parse_genes_input(DEFAULT_MVA_GENES)
        bed_path = os.path.join(out_dir, "step1_genes_bedfile.bed")
        write_bed_file(parsed_genes, bed_path)
        out_vcf_gz = os.path.join(out_dir, "step2_MVA_subset.vcf.gz")

        print(f"[Step 2] Subsetting candidate regions VCF...")
        subprocess.run(["bcftools", "view", "-R", bed_path, "-O", "z", "-o", out_vcf_gz, vcf_path], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", out_vcf_gz], check=True)
        input_vcf_gz = out_vcf_gz
    else:
        # Global genome-wide pass
        # prefer cached GTF in cache_dir to avoid re-downloads
        gtf_gz = os.path.join(cache_dir, "Homo_sapiens.GRCh38.110.gtf.gz")
        gene_bed = os.path.join(out_dir, "genes_genome_wide_pad20.bed")
        out_vcf_gz = os.path.join(out_dir, "step2_global_subset.vcf.gz")
        if not os.path.exists(gene_bed):
            if not os.path.exists(gtf_gz):
                print(f"[Step 2] Downloading Ensembl GRCh38 GTF annotation into cache dir...")
                url = "http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
                subprocess.run(["wget", "-q", "-O", gtf_gz, url], check=True)

            print(f"[Step 2] Building genic BED file from GTF (cached at {gtf_gz})...")
            intervals = []
            with gzip.open(gtf_gz, "rt") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) > 4 and parts[2] == "gene":
                        chrom = parts[0]
                        start = max(0, int(parts[3]) - 21)
                        end = int(parts[4]) + 20
                        intervals.append((chrom, start, end))

            intervals.sort(key=lambda x: (x[0], x[1], x[2]))
            merged_intervals = []
            for chrom, start, end in intervals:
                if not merged_intervals:
                    merged_intervals.append([chrom, start, end])
                else:
                    prev_c, prev_s, prev_e = merged_intervals[-1]
                    if chrom == prev_c and start <= prev_e:
                        merged_intervals[-1][2] = max(prev_e, end)
                    else:
                        merged_intervals.append([chrom, start, end])

            with open(gene_bed, "w") as f:
                for c, s, e in merged_intervals:
                    f.write(f"{c}\t{s}\t{e}\n")

        print(f"[Step 2] Subsetting PASS-only genic regions for global VEP...")
        subprocess.run(["bcftools", "view", "-f", "PASS", "-R", gene_bed, vcf_path, "-O", "z", "-o", out_vcf_gz], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", out_vcf_gz], check=True)
        input_vcf_gz = out_vcf_gz

    # Docker VEP Execution
    user_id = f"{os.getuid()}:{os.getgid()}"
    vcf_filename = os.path.basename(input_vcf_gz)

    docker_cmd = [
        "docker", "run", "--user", user_id,
        "-v", f"{cache_dir}:/opt/vep/.vep",
        "-v", f"{input_vcf_gz}:/indata/{vcf_filename}:ro",
        "-v", f"{input_vcf_gz}.tbi:/indata/{vcf_filename}.tbi:ro",
        "-v", f"{out_dir}:/out",
        "ensemblorg/ensembl-vep:release_110.1",
        "vep",
        "--input_file", f"/indata/{vcf_filename}",
        "--output_file", "/out/step2_VEP_output.txt",
        "--tab",
        "--offline", "--cache", "--dir_cache", "/opt/vep/.vep",
        "--assembly", "GRCh38", "--fork", "16", "--pick",
        "--symbol", "--show_ref_allele", "--uploaded_allele",
        "--af_gnomade", "--af_gnomadg", "--af",
        "--sift", "b", "--polyphen", "b",
        "--check_existing", "--clin_sig_allele", "1",
        "--dir_plugins", "/opt/vep/.vep/Plugins",
        "--plugin", "Phenotypes,dir=/opt/vep/.vep/Plugins/phenotype_data",
        "--gene_phenotype",
        "--fields", "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,Gene,REF_ALLELE,UPLOADED_ALLELE,SIFT,PolyPhen,AF,gnomADe_AF,gnomADg_AF,CLIN_SIG,PHENOTYPES,Existing_variation"
    ]

    print(f"[Step 2] Launching Ensembl VEP via Docker container...")
    subprocess.run(docker_cmd, check=True)
    print(f"[Step 2] VEP completed. Output saved to {vep_output}")
    return vep_output

# ----------------------------------------------------------------------
# Step 3: Rank Candidates
# ----------------------------------------------------------------------
def rank_candidates(
    vep_output_path: str,
    out_dir: str,
    method: str = "targeted",
    genes: Any = None
) -> str:
    """
    Scores and ranks candidate variants from VEP output.
    Outputs step3_ranked_candidates.csv.
    """
    os.makedirs(out_dir, exist_ok=True)
    ranked_csv = os.path.join(out_dir, "step3_ranked_candidates.csv")

    with open(vep_output_path, "r") as f:
        lines = f.readlines()

    header_idx = next(i for i, l in enumerate(lines) if l.startswith("#") and not l.startswith("##"))
    header = lines[header_idx].lstrip("#").rstrip("\n").split("\t")

    df = pd.read_csv(
        vep_output_path,
        sep="\t",
        skiprows=header_idx + 1,
        names=header,
        low_memory=False
    )

    gene_col = "SYMBOL"
    text_cols = [c for c in ["Consequence", "CLIN_SIG", "PHENOTYPES", "SIFT", "PolyPhen"] if c in df.columns]
    df[text_cols] = df[text_cols].astype("string")

    af_cols = [c for c in GNOMAD_AF_COLUMNS if c in df.columns]
    df["gnomAD_AF_max"] = df[af_cols].apply(pd.to_numeric, errors="coerce").max(axis=1).fillna(0.0)

    damaging_pattern = "|".join(CONSEQUENCE_RANK.keys())
    df["is_damaging_consequence"] = df["Consequence"].str.contains(damaging_pattern, case=False, na=False)

    rank_cols = pd.DataFrame(
        {
            term: df["Consequence"].str.contains(term, case=False, na=False) * rank
            for term, rank in CONSEQUENCE_RANK.items()
        }
    )
    df["consequence_score"] = rank_cols.max(axis=1)

    if "PHENOTYPES" in df.columns:
        df["PHENOTYPES_norm"] = normalize_phenotype_text(df["PHENOTYPES"])
        phenotype_cat_cols = pd.DataFrame({
            key: df["PHENOTYPES_norm"].str.contains(pattern, case=False, na=False, regex=True) * PHENOTYPE_WEIGHTS[key]
            for key, pattern in PHENOTYPE_PATTERNS.items()
        })
        df["phenotype_score"] = phenotype_cat_cols.sum(axis=1)
    else:
        df["phenotype_score"] = 0.0

    parsed_genes = parse_genes_input(genes)
    if method == "targeted" or parsed_genes is not None:
        candidate_set = set(parsed_genes.keys()) if parsed_genes else set(DEFAULT_MVA_GENES.keys())
        df["is_candidate_gene"] = df[gene_col].isin(candidate_set)
        df["is_rare"] = df["gnomAD_AF_max"] < RARE_AF_THRESHOLD
        filtered = df[df["is_damaging_consequence"] & df["is_rare"] & df["is_candidate_gene"]].copy()
    else:
        df["is_rare"] = df["gnomAD_AF_max"] < RARE_AF_THRESHOLD
        filtered = df[df["is_damaging_consequence"] & df["is_rare"]].copy()

    clin_col = "CLIN_SIG" if "CLIN_SIG" in df.columns else None
    filtered["has_clinvar_hit"] = filtered[clin_col].apply(has_clinvar) if clin_col else False

    sift_bad = (
        filtered["SIFT"].str.contains("deleterious", case=False, na=False)
        if "SIFT" in filtered.columns
        else False
    )
    polyphen_bad = (
        filtered["PolyPhen"].str.contains("damaging", case=False, na=False)
        if "PolyPhen" in filtered.columns
        else False
    )
    filtered["predicted_damaging"] = sift_bad | polyphen_bad

    filtered["priority_score"] = (
        filtered["has_clinvar_hit"].astype(int) * 100
        + filtered["predicted_damaging"].astype(int) * 15
        + filtered["phenotype_score"] * 10
        + filtered["consequence_score"] * 10
        + (1 - filtered["gnomAD_AF_max"]) * 5
    )

    ranked = filtered.sort_values("priority_score", ascending=False)

    display_cols = [c for c in [
        "Uploaded_variation", "Location", gene_col, "REF_ALLELE", "UPLOADED_ALLELE",
        "Consequence", "gnomAD_AF_max", "SIFT", "PolyPhen", clin_col, "PHENOTYPES", "priority_score"
    ] if c and c in ranked.columns]

    final_ranked = ranked[display_cols]

    sub_dfs = []
    for loc in set(final_ranked.Location):
        sub_df = final_ranked[final_ranked.Location == loc].copy()
        sub_df.Consequence = (",".join(set((",".join(sub_df.Consequence.to_numpy())).split(","))))
        if clin_col and clin_col in sub_df.columns:
            sub_df.CLIN_SIG = (",".join(set((",".join(sub_df.CLIN_SIG.to_numpy())).split(","))))
        if "PHENOTYPES" in sub_df.columns:
            sub_df.PHENOTYPES = (",".join(set((",".join(sub_df.PHENOTYPES.to_numpy())).split(","))))
        sub_dfs.append(sub_df)

    if sub_dfs:
        final_ranked = pd.concat(sub_dfs)
        final_ranked.sort_values("priority_score", inplace=True, ascending=False)
        final_ranked.drop_duplicates(subset=["Location", "SYMBOL"], inplace=True, keep="first")
        final_ranked.reset_index(drop=True, inplace=True)
    else:
        final_ranked = pd.DataFrame(columns=display_cols)

    final_ranked.to_csv(ranked_csv, index=False)
    print(f"[Step 3] Candidate ranking saved ({len(final_ranked)} variants) to {ranked_csv}")
    return ranked_csv

# ----------------------------------------------------------------------
# Step 4: Merge Annotations and Genotypes
# ----------------------------------------------------------------------
def merge_genotypes(
    genotypes_tsv: str,
    ranked_csv: str,
    out_dir: str
) -> str:
    """
    Merges ranked candidates with genotype calling depth and quality metrics.
    Outputs step4_merged.csv.
    """
    os.makedirs(out_dir, exist_ok=True)
    merged_csv = os.path.join(out_dir, "step4_merged.csv")

    genotypes = pd.read_csv(genotypes_tsv, sep="\t", dtype={"CHROM": str}, low_memory=False)
    candidates = pd.read_csv(ranked_csv)

    if candidates.empty:
        merged = pd.DataFrame()
        merged.to_csv(merged_csv, index=False)
        print(f"[Step 4] Candidates table empty, saved empty {merged_csv}")
        return merged_csv

    candidates[["CHROM", "POS"]] = candidates["Location"].str.extract(r"^([^:]+):(\d+)")
    candidates["POS"] = candidates["POS"].astype(int)
    genotypes["CHROM"] = genotypes["CHROM"].astype(str)

    merged = candidates.merge(
        genotypes.rename(columns={"REF": "REF_ALLELE"}),
        on=["CHROM", "POS", "REF_ALLELE"],
        how="left"
    )

    # Hard quality filters
    merged = merged[merged["FILTER"] == "PASS"]
    merged = merged[merged["GT"].isin(["0/1", "1/1", "1/0", "1|0", "0|1", "1|1"])]
    merged = merged[merged["ALT"].notna()]
    merged = merged[pd.to_numeric(merged["AD"].str.split(",").str[1], errors="coerce").fillna(0) > 0]
    merged = merged[pd.to_numeric(merged["DP"], errors="coerce") >= 10]

    merged.to_csv(merged_csv, index=False)
    print(f"[Step 4] Genotype-filtered merged results ({len(merged)} variants) saved to {merged_csv}")
    return merged_csv

# ----------------------------------------------------------------------
# Step 5: Build Hackathon Submission
# ----------------------------------------------------------------------
def _parse_variant_location(loc: str) -> Tuple[Optional[str], Optional[int]]:
    m = re.match(r"^(\w+):(\d+)(?:-(\d+))?$", str(loc))
    if not m:
        return None, None
    chrom, start, _ = m.groups()
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return chrom, int(start)

def _classify_zygosity(gt: Any) -> str:
    if pd.isna(gt) or str(gt) in (".", "./.", ".|."):
        return "no_call"
    alleles = re.split(r"[/|]", str(gt))
    if len(set(alleles)) == 1:
        return "hom_ref" if alleles[0] == "0" else "hom_alt"
    return "het"

def _build_submission_notes(row: pd.Series, extra: str = "") -> str:
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

def build_submission(
    base_dir: Optional[str] = None,
    vcf_path: str = "/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz",
    proband_id: Optional[str] = None,
    username: str = "antocad3",
    approach_name: str = "ensemble_without_alignment",
    max_rows: int = 10,
) -> str:
    """
    Builds the hackathon submission CSV from the step4_merged.csv file(s) produced
    by run()/run_all() (targeted and/or global, whichever are found under base_dir).
    Writes "{username}_{approach_name}.csv" next to this script and returns its path.
    """
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))

    step4_files = {}
    for method in ("targeted", "global"):
        path = os.path.join(base_dir, method, "step4_merged.csv")
        if os.path.exists(path):
            step4_files[method] = path
    if not step4_files:
        raise FileNotFoundError(
            f"No step4_merged.csv found under {base_dir}/targeted or {base_dir}/global. Run the pipeline first."
        )

    if proband_id is None:
        result = subprocess.run(["bcftools", "query", "-l", vcf_path], check=True, capture_output=True, text=True)
        proband_id = result.stdout.strip().splitlines()[0]

    known_mva_genes = set(DEFAULT_MVA_GENES.keys())

    frames = []
    for source_name, path in step4_files.items():
        d = pd.read_csv(path, low_memory=False)
        required = ["Location", "REF_ALLELE", "ALT", "SYMBOL", "priority_score", "GT"]
        missing = [c for c in required if c not in d.columns]
        if missing:
            raise ValueError(f"{path} is missing expected column(s): {missing}")
        d[["chrom", "pos"]] = d["Location"].apply(lambda x: pd.Series(_parse_variant_location(x)))
        d = d.dropna(subset=["chrom", "pos"]).copy()
        d["pos"] = d["pos"].astype(int)
        d["ref"] = d["REF_ALLELE"]
        d["alt"] = d["ALT"]
        d["zygosity"] = d["GT"].apply(_classify_zygosity)
        d = d[d["zygosity"].isin(["het", "hom_alt"])].copy()
        d["source"] = source_name
        print(f"[Step 5] Loaded {len(d)} candidate rows from '{source_name}' ({path})")
        frames.append(d)

    all_variants = pd.concat(frames, ignore_index=True)
    all_variants["variant_key"] = (
        all_variants["chrom"] + ":" + all_variants["pos"].astype(str) + ":"
        + all_variants["ref"].astype(str) + ">" + all_variants["alt"].astype(str)
    )

    def collapse_variant_group(grp: pd.DataFrame) -> pd.Series:
        methods = sorted(grp["source"].unique())
        best = grp.loc[grp["priority_score"].idxmax()]  # inherit annotation from best-scoring hit
        zygosity_conflict = grp["zygosity"].nunique() > 1
        consequence = ",".join(set((",".join(grp.Consequence.to_numpy())).split(",")))
        clin_sig = ",".join(set((",".join(grp.CLIN_SIG.to_numpy())).split(",")))
        phenotypes = ",".join(set((",".join(grp.PHENOTYPES.str.lower().to_numpy())).split(",")))
        return pd.Series({
            "chrom": best["chrom"], "pos": best["pos"], "ref": best["ref"], "alt": best["alt"],
            "gene": best["SYMBOL"],
            "zygosity": best["zygosity"],
            "zygosity_conflict": zygosity_conflict,
            "n_methods": len(methods),
            "methods": "+".join(methods),
            "consequence": consequence,
            "gnomad_af": best.get("gnomAD_AF_max", np.nan),
            "clin_sig": clin_sig,
            "phenotypes": phenotypes,
            "combined_score": grp["priority_score"].sum(),  # sum rewards hits found by multiple methods
        })

    variants = (
        all_variants.groupby("variant_key", group_keys=False)
        .apply(collapse_variant_group)
        .sort_values(by="combined_score", ascending=False)
        .reset_index(drop=True)
    )
    print(f"[Step 5] {len(variants)} unique variants after merging method(s) "
          f"({(variants['n_methods'] > 1).sum()} found by more than one)")

    def is_known_gene(symbol: Any) -> bool:
        return str(symbol).upper() in known_mva_genes

    # Homozygous single-variant candidates
    hom = variants[variants["zygosity"] == "hom_alt"].copy()
    hom_rows = []
    for _, r in hom.iterrows():
        hom_rows.append({
            "score": r["combined_score"],
            "chrom_1": r["chrom"], "pos_1": r["pos"], "ref_1": r["ref"], "alt_1": r["alt"],
            "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
            "is_mva_gene": is_known_gene(r["gene"]),
            "has_clinvar": has_clinvar(r["clin_sig"]),
            "notes": "gene: " + r["gene"] + " | 1st variant: " + _build_submission_notes(r),
        })

    # Compound-het candidate pairs (within each gene, top 2 by combined_score)
    het = variants[variants["zygosity"] == "het"].copy()
    het["gene"] = het["gene"].replace({"": np.nan, "-": np.nan, "None": np.nan})
    het_known = het.dropna(subset=["gene"])
    het_unknown = het[het["gene"].isna()]

    pair_rows = []
    paired_indices = set()
    for gene, grp in het_known.groupby("gene"):
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
                "gene: " + gene + " | 1st variant: " + _build_submission_notes(v1)
                + " | 2nd variant: " + _build_submission_notes(v2)
            ),
        })

    # Unpaired hets: unknown gene symbol, or the only damaging hit in their gene
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
            "notes": "gene: " + (r["gene"] if pd.notna(r["gene"]) else "<unknown>") + " | 1st variant: "
                     + _build_submission_notes(r, extra=f"single heterozygous ({reason}; "
                       f"does not independently fit the family's apparent recessive inheritance pattern)"),
        })

    candidates = pd.DataFrame(hom_rows + pair_rows)
    if candidates.empty:
        raise ValueError("No homozygous or compound-het candidates survived filtering.")

    candidates = candidates.sort_values("score", ascending=False).reset_index(drop=True)
    candidates = candidates.head(max_rows).copy()

    min_s, max_s = candidates["score"].min(), candidates["score"].max()
    if max_s > min_s:
        candidates["epcr"] = 0.05 + 0.95 * (candidates["score"] - min_s) / (max_s - min_s)
    else:
        candidates["epcr"] = 1.0

    candidates["finding_type"] = np.where(candidates["is_mva_gene"], "primary", "secondary")
    candidates["proband_id"] = proband_id

    output_cols = [
        "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
        "chrom_2", "pos_2", "ref_2", "alt_2",
        "epcr", "finding_type", "notes",
    ]
    submission = candidates[output_cols]

    output_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"{username}_{approach_name}.csv")
    submission.to_csv(output_file, index=False)
    print(f"[Step 5] Wrote {len(submission)} ranked candidate row(s) to {output_file}")
    return output_file

# ----------------------------------------------------------------------
# Top-level High-Level API Functions
# ----------------------------------------------------------------------
def run(
    method: str = "targeted",
    genes: Any = "default",
    vcf_path: str = "/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz",
    out_dir: Optional[str] = None,
    cache_dir: str = "/home/molcure/storageB/antoine/mva_hackathon/vep",
    skip_vep_if_exists: bool = False,
) -> str:
    """
    High level function to run Track 1 completely inside Python.
    
    Parameters:
        method: "targeted" or "global"
        genes: Priors for candidate genes (dict, string, list, or None). Default for targeted = 7 MVA genes. Default for global = None.
        vcf_path: Path to input VCF file
        out_dir: Directory to store step outputs. Defaults to track1/{method}
        cache_dir: Directory for VEP cache and plugins
        skip_vep_if_exists: If True, skips running docker VEP if step2_VEP_output.txt already exists.
    """
    method = str(method).strip().lower()
    if method not in ("targeted", "global"):
        raise ValueError(f"Invalid method '{method}'. Expected 'targeted' or 'global'.")

    if genes == "default":
        genes_param = DEFAULT_MVA_GENES if method == "targeted" else None
    else:
        genes_param = genes

    if out_dir is None:
        base_path = os.path.dirname(os.path.abspath(__file__))
        out_dir = os.path.join(base_path, method)

    print(f"==================================================")
    print(f"       Running Track 1 Pipeline ({method.upper()})")
    print(f"==================================================")
    print(f"VCF Path : {vcf_path}")
    print(f"Out Dir  : {out_dir}")
    print(f"Priors   : {genes_param if genes_param is not None else 'None (genome-wide)'}\n")

    # Step 1: Genotypes
    genotypes_tsv = get_genotypes(vcf_path=vcf_path, out_dir=out_dir, genes=genes_param, skip_if_exists=skip_vep_if_exists)

    # Step 2: VEP
    vep_out = run_vep(
        vcf_path=vcf_path,
        cache_dir=cache_dir,
        out_dir=out_dir,
        method=method,
        genes=genes_param,
        skip_if_exists=skip_vep_if_exists
    )

    # Step 3: Candidate Ranking
    ranked_csv = rank_candidates(
        vep_output_path=vep_out,
        out_dir=out_dir,
        method=method,
        genes=genes_param
    )

    # Step 4: Genotype Merging & Quality Filtering
    merged_csv = merge_genotypes(
        genotypes_tsv=genotypes_tsv,
        ranked_csv=ranked_csv,
        out_dir=out_dir
    )

    print(f"==================================================")
    print(f" Track 1 {method.upper()} completed successfully!")
    print(f" Result: {merged_csv}")
    print(f"==================================================\n")

    return merged_csv

def run_all(
    vcf_path: str = "/home/molcure/storageB/antoine/mva_hackathon/data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz",
    base_dir: Optional[str] = None,
    cache_dir: str = "/home/molcure/storageB/antoine/mva_hackathon/vep",
    skip_vep_if_exists: bool = False,
) -> Dict[str, str]:
    """
    Executes both targeted and global Track 1 pipelines.
    """
    if base_dir is None:
        base_dir = os.path.dirname(os.path.abspath(__file__))

    targeted_dir = os.path.join(base_dir, "targeted")
    global_dir = os.path.join(base_dir, "global")

    res_targeted = run(
        method="targeted",
        genes="default",
        vcf_path=vcf_path,
        out_dir=targeted_dir,
        cache_dir=cache_dir,
        skip_vep_if_exists=skip_vep_if_exists
    )

    res_global = run(
        method="global",
        genes=None,
        vcf_path=vcf_path,
        out_dir=global_dir,
        cache_dir=cache_dir,
        skip_vep_if_exists=skip_vep_if_exists
    )

    return {"targeted": res_targeted, "global": res_global}

# ----------------------------------------------------------------------
# CLI Interface
# ----------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Track 1 Pipeline: Targeted & Global Analysis")
    parser.add_argument("--method", choices=["targeted", "global", "all"], default="all", help="Pipeline track method")
    parser.add_argument("--vcf", default="../data/vcf/WGS_EX2312012_HGWCNDSX7.vcf.gz", help="Path to input VCF")
    parser.add_argument("--out-dir", default=None, help="Output directory")
    parser.add_argument("--cache-dir", default="/home/molcure/storageB/antoine/mva_hackathon/vep", help="VEP cache directory")
    parser.add_argument("--skip-vep-if-exists", action="store_true", help="Skip VEP docker run if output text file exists")
    parser.add_argument("--build-submission", action="store_true", help="Build the hackathon submission CSV from step4_merged.csv after running")
    parser.add_argument("--username", default="antocad3", help="Username prefix for the submission filename")
    parser.add_argument("--approach-name", default="ensemble_without_alignment", help="Approach name suffix for the submission filename")

    args = parser.parse_args()

    # Run setup once at program start to prepare cache/docker/plugins if needed
    try:
        setup_cache(args.cache_dir)
    except Exception:
        print(f"[Main] Warning: setup_cache() encountered issues; continuing anyway.")

    if args.method == "all":
        run_all(vcf_path=args.vcf, base_dir=args.out_dir, cache_dir=args.cache_dir, skip_vep_if_exists=args.skip_vep_if_exists)
    else:
        run(method=args.method, vcf_path=args.vcf, out_dir=args.out_dir, cache_dir=args.cache_dir, skip_vep_if_exists=args.skip_vep_if_exists)

    if args.build_submission:
        build_submission(base_dir=args.out_dir, vcf_path=args.vcf, username=args.username, approach_name=args.approach_name)
