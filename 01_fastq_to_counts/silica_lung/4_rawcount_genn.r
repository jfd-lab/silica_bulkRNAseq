library(Rsubread)
library(edgeR)
source("/nfs/roberts/project/pi_epm38/jfd42/pipelines/Xiting_bulkRNA_FastQtoCounts/my_functions.R")

# --- Bouchet/mouse paths ---
output.dir<-"/nfs/roberts/scratch/pi_epm38/jfd42/silica/results_Rsubread_GRCm39"
tophat.dir<-"/nfs/roberts/scratch/pi_epm38/jfd42/silica/results_STAR_GRCm39"
gtf.path<-"/nfs/roberts/project/pi_epm38/jfd42/GENCODE/GRCm39/gencode.vM39.primary_assembly.annotation.gtf"

if(file.exists(output.dir)==F){
  dir.create(output.dir)
}

dir.list<-list.files(tophat.dir)
dir.samplenames<-dir.list

# --------------------------------------------------------------------------
# FIX: original script referenced dir.samplenames[j] with no loop and no
# defined j, so it only ever (attempted to) process one undefined sample.
# All of Xiting's body variables already use [j], so we simply wrap her
# existing code in the for(j ...) loop it was written to sit inside.
# Nothing else about the featureCounts logic changes.
# --------------------------------------------------------------------------
for(j in 1:length(dir.samplenames)){
  output.filepath<-file.path(output.dir,paste0("rawcount_",dir.samplenames[j],".txt"))

  sample.name<-dir.samplenames[j]
  bam.filepath<-list.files(file.path(tophat.dir,dir.samplenames[j]))
  bam.filepath<-bam.filepath[grep("Aligned.sortedByCoord.out.bam",bam.filepath)]
  bam.filepath<-file.path(tophat.dir,dir.samplenames[j],bam.filepath)

  # Single-end: isPairedEnd=F, requireBothEndsMapped=F
  temp<-featureCounts(bam.filepath,annot.ext=gtf.path,isGTFAnnotationFile=T,GTF.featureType="exon",GTF.attrType="gene_name",isPairedEnd=F,useMetaFeatures=T,allowMultiOverlap=T,nthreads=1,countMultiMappingReads=T)
  cmd.out<-cbind(temp$annotation,temp$counts)

  write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=T,quote=F)
}
