#!/bin/bash
# =============================================================================
# 9_setup_sh_matching.sh
# One-time setup for local UNITE Species Hypothesis (SH) matching, used by
# 9_dark_taxa_SH_matching.R for the E1 "dark taxa" exploration.
#
# Builds the authoritative TU-NHM sh_matching pipeline (the same service behind
# PlutoF / unite.ut.ee) as a Singularity/Apptainer container and prepares the
# UNITE SH reference UDB databases, so SH assignment can be run reproducibly on
# the server rather than via the manual web upload.
#
# Reference: github.com/TU-NHM/sh_matching_pub, pinned to release v2.0.4.
# Cite: Abarenkov K, Koljalg U, Nilsson RH (2022) UNITE Species Hypotheses
#       Matching Analysis. Biodiv. Inf. Sci. Stand. 6: e93856.
#
# Idempotent: safe to re-run; steps already completed are skipped.
# Run once:  bash Scripts/9_setup_sh_matching.sh
# =============================================================================
set -euo pipefail

BASE=/home/daniel/Ptarmigan/sh_matching
REPO="$BASE/sh_matching_pub"
SIF="$BASE/sh_matching.sif"
TAG=v2.0.4
# UNITE SH UDB reference package (PlutoF file repository view/8135714)
UDB_URL="https://s3.hpc.ut.ee/plutof-public/original/b4e6594e-c4bc-42b7-a4c6-8286e8b73943.zip"

mkdir -p "$BASE"
cd "$BASE"

# ---- 1. Get pipeline source, pinned to a release tag (master is dev) --------
if [ ! -d "$REPO" ]; then
  echo "[1/5] Cloning sh_matching_pub ..."
  git clone --depth 1 --branch "$TAG" https://github.com/TU-NHM/sh_matching_pub.git
else
  echo "[1/5] Repo present; ensuring tag $TAG ..."
  ( cd "$REPO" && git fetch --depth 1 origin tag "$TAG" && git checkout -q "$TAG" )
fi

# ---- 2. Working directories the pipeline expects in the run cwd -------------
echo "[2/5] Creating indata/ outdata/ userdir/ data_udb/ ..."
mkdir -p indata outdata userdir data_udb

# ---- 3. Build the container (needs user-namespace fakeroot; verified OK on
#         this host: unprivileged_userns_clone=1 + subuid/subgid for daniel) --
if [ ! -f "$SIF" ]; then
  echo "[3/5] Building $SIF (this downloads USEARCH/VSEARCH/ITSx/Krona + an in-image UNITE db; slow) ..."
  ( cd "$REPO" && singularity build --fakeroot "$SIF" sh_matching.def )
else
  echo "[3/5] $SIF already built; skipping."
fi

# ---- 4. Download + unpack the SH reference FASTAs (~144 MB zip) -------------
if [ ! -f data_udb/sanger_refs_sh.fasta ] && [ ! -f data_udb/sanger_refs_sh.udb ]; then
  echo "[4/5] Downloading UNITE SH reference package ..."
  wget -q "$UDB_URL" -O sh_matching_data_udb_0_5.zip
  unzip -o -q sh_matching_data_udb_0_5.zip
  rm -f sh_matching_data_udb_0_5.zip
else
  echo "[4/5] Reference FASTAs/UDBs already present; skipping download."
fi

# ---- 5. Build the UDB indexes. The host has no vsearch, so use the copy
#         bundled inside the container. ---------------------------------------
VSEARCH=/sh_matching/programs/vsearch/bin/vsearch
if [ ! -f data_udb/sanger_refs_sh.udb ]; then
  echo "[5/5] Building sanger_refs_sh.udb ..."
  singularity exec "$SIF" "$VSEARCH" --makeudb_usearch data_udb/sanger_refs_sh.fasta      --output data_udb/sanger_refs_sh.udb
fi
if [ ! -f data_udb/sanger_refs_sh_full.udb ]; then
  echo "[5/5] Building sanger_refs_sh_full.udb ..."
  singularity exec "$SIF" "$VSEARCH" --makeudb_usearch data_udb/sanger_refs_sh_full.fasta --output data_udb/sanger_refs_sh_full.udb
fi

echo
echo "Setup complete. Artifacts:"
ls -la "$SIF" data_udb/*.udb 2>/dev/null || true
echo "Next: conda run -n r_env Rscript Scripts/9_dark_taxa_SH_matching.R"
