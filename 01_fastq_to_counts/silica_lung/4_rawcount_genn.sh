#!/bin/bash
# ============================================================================
# 4_rawcount_genn.sh
# SLURM wrapper to run 4_rawcount_genn.r (featureCounts, Rsubread) on Bouchet.
#
# featureCounts over ~30 BAMs against the full mouse GTF is real work and must
# NOT run on a login node. Submit this with:  sbatch 4_rawcount_genn.sh
#
# Uses the Bioconductor bundle, which provides R + Rsubread + edgeR together
# (matches Bouchet's 2024a toolchain).
# ============================================================================
#SBATCH --partition=day
#SBATCH --job-name=rawcount_GRCm39
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=NONE
#SBATCH --error=rawcount_GRCm39.e%J
#SBATCH --output=rawcount_GRCm39.o%J

module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2

# Path to the R script (lives in the pipeline code folder alongside the other scripts)
Rscript /nfs/roberts/project/pi_epm38/jfd42/pipelines/Xiting_bulkRNA_FastQtoCounts/4_rawcount_genn.r
