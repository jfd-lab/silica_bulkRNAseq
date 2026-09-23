# Silica bulk RNA-seq pipeline — Bouchet / mouse port

Ported from Xiting Yan's `2026_Manning_PAvRVBulk` pipeline (human GRCh38, McCleary HPC)
to your silica Plasmidsaurus SS8VZG data (mouse GRCm39, Bouchet HPC).

## What changed vs. Xiting's original (and ONLY these things)

1. **Human → mouse.** GRCh38/GENCODE v50 → GRCm39/GENCODE M39. New reference download
   script (`00_download_reference.sh`). All `--genomeDir`, `--sjdbGTFfile`, and
   featureCounts GTF paths point at the mouse M39 primary-assembly FASTA + GTF.
2. **McCleary → Bouchet.** All `/home/xy48/scratch_pi_xy48/xy48/...` paths → your Roberts
   paths (`/nfs/roberts/{home,project,scratch}/pi_epm38/jfd42/...`). Partition
   `day,week,day_amd` → `day` (McCleary's `day_amd` does not exist on Bouchet). Module
   strings confirmed present on Bouchet and unchanged (STAR 2.7.11b, cutadapt 5.1,
   FastQC 0.12.1; Subread 2.1.1 for step 4).
3. **rawcount loop bug fixed.** The original `rawcount_genn.r` referenced
   `dir.samplenames[j]` with no loop and undefined `j`. Wrapped the existing body in the
   `for(j in 1:length(dir.samplenames))` loop it was written for. **No other logic changed** —
   `gene_name` counting, `allowMultiOverlap=T`, `countMultiMappingReads=T` all kept as-is.

Everything else (STAR flags, 3-round Cutadapt, featureCounts settings) is byte-for-byte Xiting.

## Folder layout (code vs. data)

Two separate locations, on purpose:

**Code folder** (persistent project space) — put ALL of these together here:
```
/nfs/roberts/project/pi_epm38/jfd42/pipelines/Xiting_bulkRNA_FastQtoCounts/
    00_download_reference.sh
    1_FastQC.Rmd
    2_CutAdapt.Rmd
    3_STARmapping.Rmd
    4_rawcount_genn.r
    4_rawcount_genn.sh
    5_merge_counts.r
    my_functions.R      <-- copy Xiting's file in here yourself
    README.md
```
Every `source(...)` line points at `my_functions.R` in this folder, so copying it in is all that's needed.

**Data / outputs** (referenced BY the scripts, not stored in the code folder):
- fastqs (input):   `/nfs/roberts/project/pi_epm38/jfd42/silica/plasmidsaurus_SS8VZG_raw_data/FASTQ_files/{lung,rv}/`
- reference genome:  `/nfs/roberts/project/pi_epm38/jfd42/GENCODE/GRCm39/`  and STAR index in `.../STARindex/GRCm39_94/`
- working outputs:   `/nfs/roberts/scratch/pi_epm38/jfd42/silica/`  (scratch — big BAMs/intermediates live here; purged after 60 days)

Keep big data OUT of the code folder — that's why `work.dir` points at scratch, not here.

## Run order

0. `bash 00_download_reference.sh`  — download GRCm39 FASTA + GTF to project space (once).
   Then build the STAR index (command block at top of `3_STARmapping.Rmd`).
1. `1_FastQC.Rmd`      — raw-read QC. Generates per-sample `.sh`, run `jobsub.bat`.
2. `2_CutAdapt.Rmd`    — 3-round trim (Illumina + polyA + polyG) + post-trim FastQC.
3. `3_STARmapping.Rmd` — STAR two-pass mapping + mapping summary table.
4. `4_rawcount_genn.r` (submitted via `4_rawcount_genn.sh`, i.e. `sbatch 4_rawcount_genn.sh`) —
   featureCounts per sample → `results_Rsubread_GRCm39/rawcount_<sample>.txt`. Do NOT run on a
   login node; it's real compute.
5. `5_merge_counts.r` — merge per-sample files → single `rawcounts_GRCm39.txt` (annotation cols
   1-6 + one count column per sample), the exact format the downstream DESeq2 code reads. This
   is what drops into your Jupyter QC/DESeq2 notebooks.

**Lung and RV are separate runs.** Each script has a `data.dir` pointing at the lung folder;
for the RV analysis, change `data.dir` to your RV fastq folder and rerun. Outputs are named
generically, so either keep separate `work.dir`s per tissue or rename outputs between runs so
lung/RV don't overwrite each other.

## MUST-VERIFY before / on first run

- [ ] **Read length / sjdbOverhang.** Scripts default to Xiting's `sjdbOverhang=93` (= 94 bp − 1).
      Confirm your Plasmidsaurus read length:
      `zcat <one_sample>.fastq.gz | head -400 | awk 'NR%4==2{print length($0)}' | sort -n | tail -1`
      If not 94 bp, change `sjdbOverhang` in BOTH the index-build block and the STAR call in
      `3_STARmapping.Rmd`, and rename the index dir accordingly.
- [ ] **R + Rsubread/edgeR for step 4 (CONFIRMED available).** Bouchet has
      `R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2`, which provides R + Rsubread + edgeR
      together. `4_rawcount_genn.sh` loads it. Quick sanity check that the two packages are in
      the bundle:
      `module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 && Rscript -e 'library(Rsubread);library(edgeR);cat("both OK\n")'`
      If either is missing, `BiocManager::install(...)` it into a personal library. Steps 1–3
      need no R on the compute node.
- [ ] **my_functions.R location.** All scripts source
      `/nfs/roberts/project/pi_epm38/jfd42/pipelines/Xiting_bulkRNA_FastQtoCounts/my_functions.R`.
      Copy Xiting's `my_functions.R` into that folder (alongside these scripts) before running.
      It provides `my.element.extract`, etc. Nothing else references her old repo.
- [ ] **fastq folder path.** Confirm the lung/RV fastq paths under
      `/nfs/roberts/project/pi_epm38/jfd42/silica/plasmidsaurus_SS8VZG_raw_data/FASTQ_files/`.
      Note: `/home/jfd42/project_pi_epm38/...` is likely a symlink to `/nfs/roberts/project/...`;
      use `readlink -f` and prefer the canonical `/nfs/roberts/...` path inside SLURM jobs.
- [ ] **Sample IDs vs. metadata.** Filename-before-first-dot becomes the sample ID. Confirm
      these match the keys your Jupyter QC/DESeq2 notebooks + metadata.csv expect (watch the
      FC4/4FC ambiguity you flagged).

## Merging per-sample counts → one matrix

Handled by `5_merge_counts.r` (now scripted, no manual step). It takes the annotation columns
once and cbinds each sample's count column with the sample name as header, producing
`rawcounts_GRCm39.txt` — annotation in cols 1-6, one count column per sample, exactly the layout
Xiting's `4_DESeq2.Rmd` reads (`annot <- [,1:6]; counts <- [,7:ncol]`). It also asserts gene
order matches across files before merging. Drop the result into your existing Jupyter QC
(`load_plasmidsaurus`) and DESeq2 notebooks and re-run the whole analysis on the new counts (no
need to compare against the old Plasmidsaurus matrix).

## Note on strandedness (flag only — nothing changed)

featureCounts here runs unstranded (Xiting's setting; no `strandSpecific` argument).
Plasmidsaurus 3'-tag libraries are usually stranded. This is kept as-is for methodological
consistency and won't break anything. If the featureCounts assignment rate looks low in the
per-sample summary, testing `strandSpecific=1` (or `=2`) on one sample is the first thing to check.
