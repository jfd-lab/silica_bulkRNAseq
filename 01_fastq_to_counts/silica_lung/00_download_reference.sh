#!/bin/bash
# ============================================================================
# 00_download_reference.sh
# Download the GENCODE mouse GRCm39 (Release M39) reference genome + annotation
# for the silica bulk RNA-seq pipeline on the Bouchet HPC.
#
# This is the mouse analogue of the human reference Xiting used:
#   human FASTA : GRCh38.primary_assembly.genome.fa
#   human GTF   : gencode.v50.primary_assembly.annotation.gtf
#   mouse FASTA : GRCm39.primary_assembly.genome.fa        (this script)
#   mouse GTF   : gencode.vM39.primary_assembly.annotation.gtf (this script)
#
# Run this ONCE on a Bouchet login node (or an interactive devel session).
# The reference lives on PROJECT space (persistent; scratch is purged after 60 days).
#
# Usage:
#   bash 00_download_reference.sh
# ============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------
# Destination: persistent project space (NOT scratch).
# From `mydirectories`:  project: /nfs/roberts/project/pi_epm38/jfd42
# ----------------------------------------------------------------------------
REF_DIR="/nfs/roberts/project/pi_epm38/jfd42/GENCODE/GRCm39"
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "Downloading GENCODE mouse Release M39 (GRCm39) into: ${REF_DIR}"

# ----------------------------------------------------------------------------
# 1. Primary-assembly genome FASTA (chromosomes + scaffolds)
#    ~800 MB gzipped. Primary assembly is the correct input for alignment.
# ----------------------------------------------------------------------------
wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/GRCm39.primary_assembly.genome.fa.gz"

# ----------------------------------------------------------------------------
# 2. Primary-assembly comprehensive annotation GTF
# ----------------------------------------------------------------------------
wget -c "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39/gencode.vM39.primary_assembly.annotation.gtf.gz"

# ----------------------------------------------------------------------------
# Decompress (STAR needs the plain .fa and .gtf, not gzipped)
# ----------------------------------------------------------------------------
echo "Decompressing..."
gunzip -k GRCm39.primary_assembly.genome.fa.gz
gunzip -k gencode.vM39.primary_assembly.annotation.gtf.gz

echo ""
echo "Done. Reference files:"
echo "  FASTA : ${REF_DIR}/GRCm39.primary_assembly.genome.fa"
echo "  GTF   : ${REF_DIR}/gencode.vM39.primary_assembly.annotation.gtf"
echo ""
echo "Next: build the STAR index (see 3_STARmapping.Rmd, 'Build Genome Index' section)."
