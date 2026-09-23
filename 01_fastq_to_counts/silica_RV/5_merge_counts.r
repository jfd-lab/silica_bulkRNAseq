# ============================================================================
# 5_merge_counts.r
# Merge the per-sample rawcount_<sample>.txt files produced by
# 4_rawcount_genn.r into ONE count matrix:
#   columns 1-6 = featureCounts annotation (GeneID, Chr, Start, End, Strand, Length)
#   columns 7+  = one integer count column per sample (header = sample name)
#
# This is the format Xiting's downstream DESeq2 code reads:
#   annot.matrix <- count.matrix[,1:6]; count.matrix <- count.matrix[,7:ncol]
#
# Run after 4_rawcount_genn.r finishes. Not compute-heavy; a login node is fine,
# but you can wrap it like 4_rawcount_genn.sh if you prefer.
#
#   module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
#   Rscript 5_merge_counts.r
# ============================================================================

count.dir<-"/nfs/roberts/scratch/pi_epm38/jfd42/silica_RV/results_Rsubread_GRCm39"
output.filepath<-"/nfs/roberts/scratch/pi_epm38/jfd42/silica_RV/rawcounts_GRCm39_RV.txt"

files<-list.files(count.dir, pattern="^rawcount_.*\\.txt$", full.names=FALSE)
if(length(files)==0){
  stop("No rawcount_*.txt files found in ", count.dir)
}

# Sample name = filename minus the "rawcount_" prefix and ".txt" suffix
sample.names<-sub("^rawcount_","",sub("\\.txt$","",files))

merged<-NULL
annot<-NULL

for(i in 1:length(files)){
  tab<-read.table(file.path(count.dir,files[i]),sep="\t",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
  # featureCounts single-sample output: 6 annotation cols + 1 count col (col 7)
  this.annot<-tab[,1:6]
  this.count<-tab[,7]

  if(is.null(merged)){
    annot<-this.annot
    merged<-data.frame(this.count, check.names=FALSE)
    colnames(merged)<-sample.names[i]
  } else {
    # Safety: confirm gene order matches the first file before cbind
    if(!identical(annot[,1], this.annot[,1])){
      stop("Gene order mismatch in ", files[i], " vs first file. Reorder before merging.")
    }
    merged[[sample.names[i]]]<-this.count
  }
}

out<-cbind(annot, merged)
write.table(out, file=output.filepath, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

cat("Merged", length(files), "samples ->", output.filepath, "\n")
cat("Dimensions:", nrow(out), "genes x", length(files), "samples\n")
