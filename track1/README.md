# Track 1 Pipeline (MVA Hackathon)

Genotype extraction, dockerized Ensembl VEP annotation, candidate ranking, and
genotype merging for the MVA candidate gene analysis.

> [Submitted Report](antocad3_track1_report.md)
>
> [Full Method Description](antocad3_methods_description.md)


## Prerequisites

Install these before running `track1.py`:

- **Docker** (engine installed, daemon running, and your user able to run
  `docker` without `sudo` — e.g. member of the `docker` group)
- **bcftools** and **tabix** (from htslib/samtools) — used directly on the
  host to subset and index VCFs
- **git** — used to clone the `VEP_plugins` repo into the cache directory
- **wget** — only needed for `--method global`, to fetch the Ensembl GTF
- **Python 3** with the packages in `requirements.txt`:
  ```
  pip install -r requirements.txt
  ```

On Ubuntu/Debian, the system tools can be installed with:
```
sudo apt-get update
sudo apt-get install -y docker.io bcftools tabix git wget
sudo usermod -aG docker "$USER"   # then log out/in (or `newgrp docker`)
```

## First run

The first run of `python3 track1.py` will automatically:
- pull the `ensemblorg/ensembl-vep` Docker image,
- download the VEP offline cache (~20GB, can take a long time),
- clone the `VEP_plugins` repo,
- pre-fetch the Phenotypes plugin data,

into the cache directory (default: `/home/molcure/storageB/antoine/mva_hackathon/vep`).
Subsequent runs skip all of the above if already present.

## Usage

```
python3 track1.py                          # run both targeted and global methods
python3 track1.py --method targeted        # run one method only
python3 track1.py --skip-vep-if-exists      # skip re-running steps 1-2 if their outputs already exist
python3 track1.py --cache-dir /path/to/vep  # use a different VEP cache directory
```
