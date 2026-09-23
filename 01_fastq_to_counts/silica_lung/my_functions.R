my.104.colors<-c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
"#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
"#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
"#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
"#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
"#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
"#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
"#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
"#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
"#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
"#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
"#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
"#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")
#-----------------------------------------------------------------------------------------------------
# this function loads the lines in a matrix file into a matrix. First line must be the column names and there should be no row names.
#-----------------------------------------------------------------------------------------------------

my.matrix.load<-function(x){
# this function load in the motif occurence.bed file into a matrix
temp<-readLines(x)
temp.names<-unlist(strsplit(temp[1],split="\t"))
temp<-temp[-1]
temp<-strsplit(temp,split="\t")

for(i in 1:length(temp)){
if(length(temp[[i]])<length(temp.names)){
temp[[i]]<-c(temp[[i]],rep("",length(temp.names)-length(temp[[i]])))
}
}
temp<-matrix(unlist(temp),nrow=length(temp),byrow=T)
colnames(temp)<-temp.names
return(temp)
}


#-----------------------------------------------------------------------------------------------------
# this function loads the lines in a bed file into a matrix
#-----------------------------------------------------------------------------------------------------

my.bed.load<-function(x){
temp<-readLines(x)
temp<-temp[-1]
temp<-matrix(unlist(strsplit(temp,split="\t")),nrow=length(temp),byrow=T)
return(temp)
}


#-----------------------------------------------------------------------------------------------------
# this function count the original read number, unmapped read number and the number of hits in BAM files by tophat
my.tophat.count<-function(fastq.filepath,tophat.dir,output.filepath,readnum.perpair=1){

temp<-system(paste("wc -l ",fastq.filepath,sep=""),intern=TRUE)
if(readnum.perpair==2){
fastq.num<-as.numeric(my.element.extract(temp,splitchar=" ",index=1))/2
}else{
fastq.num<-as.numeric(my.element.extract(temp,splitchar=" ",index=1))/4
}
unmapped.filepath<-file.path(tophat.dir,"unmapped.bam")

temp<-system(paste("samtools flagstat ",unmapped.filepath,sep=""),intern=TRUE)
temp1<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=1))
temp2<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=3))
unmapped.num<-temp1+temp2

mapped.filepath<-file.path(tophat.dir,"accepted_hits.bam")
temp<-system(paste("samtools flagstat ",mapped.filepath,sep=""),intern=TRUE)
temp1<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=1))
temp2<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=3))
bam.num<-temp1+temp2

result<-c(fastq.num,unmapped.num,bam.num)
cat(result,file=output.filepath,append=F,sep="\n")
return(0)
}



#-----------------------------------------------------------------------------------------------------
# this function count the original read number, unmapped read number and the number of hits in BAM files by bowtie2
my.bowtie2.count<-function(fastq.filepath,bam.filepath,output.filepath,readnum.perpair=1){
temp<-system(paste("wc -l ",fastq.filepath,sep=""),intern=TRUE)
if(readnum.perpair==2){
fastq.num<-as.numeric(my.element.extract(temp,splitchar=" ",index=1))/2
}else{
fastq.num<-as.numeric(my.element.extract(temp,splitchar=" ",index=1))/4
}
#unmapped.filepath<-file.path(tophat.dir,"unmapped.bam")

temp<-system(paste("samtools flagstat ",bam.filepath,sep=""),intern=TRUE)
temp1<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=1))
temp2<-as.numeric(my.element.extract(temp[1],splitchar=" ",index=3))
total.num<-temp1+temp2

temp1<-as.numeric(my.element.extract(temp[7],splitchar=" ",index=1))
temp2<-as.numeric(my.element.extract(temp[7],splitchar=" ",index=3))
mapped.num<-temp1+temp2

bam.num<-as.numeric(my.element.extract(temp[3],splitchar=" ",index=1))
bam.num<-bam.num+as.numeric(my.element.extract(temp[3],splitchar=" ",index=3))

result<-c(fastq.num,mapped.num,mapped.num/total.num,bam.num)
cmd.out<-cbind(c("FASTQ_num","MAPPED_num","MAPPED_rate","BAM_num"),result)
write.table(cmd.out,file=output.filepath,append=F,sep="\t",row.names=F,col.names=F,quote=F)
return(0)
}

#-----------------------------------------------------------------------------------------------------
# this function replaces the every character of orig with replace
#-----------------------------------------------------------------------------------------------------
my.char.replace<-function(x,orig="_",replace="-"){
result<-unlist(strsplit(x,split=orig))
result<-paste(result,collapse=replace)
return(unname(result))
}
#-----------------------------------------------------------------------------------------------------
my.list.element.extract<-function(x,index=1){
if(index>0){
return(x[index])
}else{
return(x[length(x)])
}
}
#-----------------------------------------------------------------------------------------------------
# this function split x by splitchar and return the parts with index=index
my.element.extract<-function(x,splitchar="\t",index=1){
temp<-unlist(strsplit(x,split=splitchar))
if(index>length(temp)){
cat("there are not enough elements!\n")
return(NA)
}
if(index<0){
result<-temp[length(temp)]
}else{
result<-temp[index]
}
return(unname(result))
}


#-----------------------------------------------------------------------------------------------------
# This function remote the element at position index in x seperated by splitchar
my.element.remove<-function(x,splitchar="",index=-1){
temp<-unlist(strsplit(x,split=splitchar))
if(index>length(temp)){
cat("there are not enough elements!\n")
return(NA)
}

if(index<0){
temp<-temp[-length(temp)]
}else{
temp<-temp[-index]
}
if(splitchar=="\\."){
return(unname(paste(temp,collapse=".")))
}else{
return(unname(paste(temp,collapse=splitchar)))
}
}



#-----------------------------------------------------------------------------------------------------
# this function generate the array id from the file names of the cel file
my.celfile2arrayid<-function(x){

temp<-unlist(strsplit(x,split="\\."))[1]
temp<-unlist(strsplit(temp,split="_"))
temp<-toupper(temp)

if(substr(temp[3],nchar(temp[3])-1,nchar(temp[3]))%in%c("-S","-B")){
temp[3]<-substr(temp[3],1,nchar(temp[3])-2)
}

#temp<-paste(temp[3],"_",temp[4],sep="")
return(c(temp[3],temp[4]))
}


#-----------------------------------------------------------------------------------------------------
# This function loads in the clinical data and remove the visits with no arrays hybridized
my.clinic.dataload<-function(x,perl.path=""){
	# this function loads in the clinical data
	library(gdata)
	if(perl.path==""){
		clinic.data<-as.matrix(read.xls(x,check.names=F))
	}else{
		clinic.data<-as.matrix(read.xls(x,check.names=F,perl=perl.path))
	}
	clinic.data<-clinic.data[,colnames(clinic.data)!=""]
	cat("There are ",sum(clinic.data[,"YMD Sample Name - B"]=="" & clinic.data[,"YMD Sample Name - S"]==""),"visits with no arrays!\n")
	# remove visits that do not have arrays
	clinic.data<-clinic.data[clinic.data[,"YMD Sample Name - B"]!="" | clinic.data[,"YMD Sample Name - S"]!="",]
	return(clinic.data)
}

#-----------------------------------------------------------------------------------------------------
# this function change the date format from "mmddyy" to "yymmdd"
my.date.format<-function(x){
temp<-unlist(strsplit(x,split=""))
result<-temp[5:6]
result<-c(result,temp[1:4])
result<-paste(result,collapse="")
return(result)
}





#-----------------------------------------------------------------------------------------------------
# this function loads in the data generated by XPS R package. 
# The sample names have the format of "[B/S for tissue]_[arrayID]_[scandate].CEL"
# no detection p-values are included;
#-----------------------------------------------------------------------------------------------------
my.xps.rma.load<-function(x){
my.data<-as.matrix(read.table(x,header=T,row.names=1,sep="\t"))
if(colnames(my.data)[1]=="UNIT_ID"){
my.data<-my.data[,-1]
}
sample.names<-colnames(my.data)
scan.date<-unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1))
scan.date<-unname(sapply(scan.date,my.element.extract,splitchar="_",index=-1))
#scan.date<-sapply(scan.date,my.date.format)
array.id<-toupper(unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1)))
array.id<-unname(sapply(array.id,my.element.remove,splitchar="_",index=1))
array.id<-unname(sapply(array.id,my.element.remove,splitchar="_",index=-1))
array.id<-unname(sapply(array.id,my.char.replace,orig="_",replace="-"))
# generate the data matrix
data.matrix<-my.data
colnames(data.matrix)<-array.id
result<-list(DataMatrix=data.matrix,ScanDate=scan.date)
return(result)
}





#-----------------------------------------------------------------------------------------------------
# this function loads in the data generated by XPS R package. 
# The sample names have the format of "[B/S for tissue]_[arrayID]_[scandate].CEL"
# no detection p-values are included;
#-----------------------------------------------------------------------------------------------------
my.xps.dabg.load<-function(x){
my.data<-as.matrix(read.table(x,header=T,sep="\t"))
#if(colnames(my.data)[1]=="ProbeID"| colnames(my.data)[1]=="UNIT_ID"){
my.data<-my.data[,-1]
rownames(my.data)<-my.data[,1]
my.data<-my.data[,-1]
#}

sample.names<-colnames(my.data)
scan.date<-unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1))
scan.date<-unname(sapply(scan.date,my.element.extract,splitchar="_",index=-1))
#scan.date<-sapply(scan.date,my.date.format)
array.id<-toupper(unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1)))
array.id<-unname(sapply(array.id,my.element.remove,splitchar="_",index=1))
array.id<-unname(sapply(array.id,my.element.remove,splitchar="_",index=-1))
array.id<-unname(sapply(array.id,my.char.replace,orig="_",replace="-"))
# generate the data matrix
data.matrix<-my.data
colnames(data.matrix)<-array.id
result<-list(DataMatrix=data.matrix,ScanDate=scan.date)
return(result)
}






#-----------------------------------------------------------------------------------------------------
# this function loads in the detection p-value matrix generated by XPS R package. 
# The sample names have the format of "[B/S for tissue]_[arrayID]_[scandate].CEL"
# no detection p-values are included;
#-----------------------------------------------------------------------------------------------------
my.dabg.load<-function(x){
my.data<-as.matrix(read.table(x,header=T,sep="\t"))
# get rid of the first column "UNIT_ID"
my.data<-my.data[,-1]
rownames(my.data)<-my.data[,1]
# get rid of the second column "Probe ID"
my.data<-my.data[,-1]
# generate the names for all the samples
sample.names<-colnames(my.data)

# generate the Scan Dates for the samples
scan.date<-unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1))
scan.date<-unname(sapply(scan.date,my.element.extract,splitchar="_",index=-1))
scan.date<-sapply(scan.date,my.date.format)

# generate the array IDs from the sample names
temp<-toupper(unname(sapply(sample.names,my.element.extract,splitchar="\\.",index=1)))
temp<-toupper(unname(sapply(temp,my.element.remove,splitchar="_",index=1)))
temp<-toupper(unname(sapply(temp,my.element.remove,splitchar="_",index=-1)))
array.id<-temp
array.id<-unname(sapply(array.id,my.char.replace,orig="_",replace="-"))
# generate the data matrix
colnames(data.matrix)<-array.id
result<-list(DataMatrix=data.matrix,ScanDate=scan.date)
return(result)
}



#-----------------------------------------------------------------------------------------------------
# this function loads in the gene expression data matrix adjusted by ComBat
#-----------------------------------------------------------------------------------------------------
my.combat.load<-function(x){
my.data<-as.matrix(read.table(x,header=T,sep="\t",row.names=1,check.names=F))
if(colnames(my.data)[1]=="ProbeID"){
my.data<-my.data[,-1]
}
return(my.data)
}


#-----------------------------------------------------------------------------------------------------
# this function extract the list of entrez ID from the given gene_assignment in the annotation file

my.cor.test<-function(x,y,method.name="spearman"){
	return(list(pvalue=cor.test(x,y,alternative="two.sided",method=method.name,na.action="na.omit")$p.value,cor=cor.test(x,y,alternative="two.sided",method=method.name,na.action="na.omit")$estimate))
}



#-----------------------------------------------------------------------------------------------------
# this function extract the list of entrez ID from the given gene_assignment in the annotation file

my.entrez.extract<-function(x){
	result<-character()
	temp<-unlist(strsplit(x,split=" /// "))
	for(i in 1:length(temp)){
		temp1<-unlist(strsplit(temp[i],split=" // "))
		result<-c(result,temp1[length(temp1)])		
	}
	result<-unique(result)
	return(result)
}

#-----------------------------------------------------------------------------------------------------
# this function extract the list of entrez ID from the given gene_assignment in the annotation file
my.genesymbol.extract<-function(x){
	result<-character()
	temp<-unlist(strsplit(x,split=" /// "))
	for(i in 1:length(temp)){
		temp1<-unlist(strsplit(temp[i],split=" // "))
		result<-c(result,temp1[2])		
	}
	result<-unique(result)
	if(length(result)==1){
		return(result)
	}else{
		return(paste(result,collapse=";"))
	}
}
#-----------------------------------------------------------------------------------------------------
# this function extract the list of entrez ID from the given gene_assignment in the annotation file
my.geneid.extract<-function(x){
	result<-character()
	temp<-unlist(strsplit(x,split=" /// "))
	for(i in 1:length(temp)){
		temp1<-unlist(strsplit(temp[i],split=" // "))
		#result<-c(result,temp1[length(temp1)])	
		result<-c(result,temp1[2])		
	}
	result<-unique(result)
	return(result)
}

#-----------------------------------------------------------------------------------------------------
# this function calculate the p-value by fisher's exact test for gene set enrichemtn
my.fisherGSEA<-function(genesetgenes,selectgenes,universegenes){
temp1<-rep(0,length(universegenes))
temp2<-rep(0,length(universegenes))
names(temp1)<-universegenes
names(temp2)<-universegenes
temp1[selectgenes]<-1
temp2[genesetgenes]<-1
cat(pathwayname,"\n",sep="")
return(fisher.test(as.matrix(table(temp1,temp2)))$p.value)
}


#-----------------------------------------------------------------------------------------------------
# this function construct a vector rep(x[1],as.numeric(x[2]))
my.vect.construct<-function(x){
return(unname(rep(x[1],as.numeric(x[2]))))
}
#-----------------------------------------------------------------------------------------------------
# this function generate a list where each element contains the probeset IDs for the given pathway
my.list.construct<-function(x,entrezid2probeset){
	temp<-unique(unname(entrezid2probeset[entrezid2probeset[,1]%in%x,2]))
	return(temp)
}

#-----------------------------------------------------------------------------------------------------
# this function replace the characters "()/" with "" in the given x
my.pathway.namestand<-function(x,replace="_"){
temp<-unlist(strsplit(x,split=""))
temp[temp=="("]<-replace
temp[temp==")"]<-replace
temp[temp=="/"]<-replace
return(paste(temp,collapse="",sep=""))
}

#-----------------------------------------------------------------------------------------------------
# this function loads in gene lists pathways in MsigDB and the annotation of affymetrix probesets
# The function will return a list, with each element corresponding to each pathway and containing
# the probe sets on the array for this pathways in MsigDB. We may change the name of the pathways 
# a little bit so that the names can be used as a file name
#-----------------------------------------------------------------------------------------------------
my.pathwaylist.extract<-function(msigdb.filepath,anno.filepath){
##################
# load in the gene list for each MsigDB pathway
temp<-readLines(msigdb.filepath)
temp.names<-unname(sapply(temp,my.element.extract,splitchar="\t",index=1))
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp<-unname(sapply(temp,my.element.remove,splitchar="\t",index=1))
temp<-strsplit(temp,split="\t")
names(temp)<-temp.names
pathway2entrezid<-temp

##################
# load in the mapping from probe set names to entrez ids
#anno.filepath<-"/home/xy48/Asthma/Data/GeneExpression/ArrayAnnotation/HuGene-1_0-st-v1.na31.hg19.transcript.csv"
array.anno<-as.matrix(read.csv(anno.filepath,check.names=F,comment.char="#"))
#temp<-array.anno[array.anno[,8]!="---",]
temp<-array.anno[array.anno[,"category"]=="main" & array.anno[,"gene_assignment"]!="---",]
temp.entrez<-sapply(temp[,8],my.entrez.extract)
names(temp.entrez)<-temp[,1]
probeset2entrezid<-temp.entrez

####################
# generate the mapping from entrez ID to probe set ID
temp.length<-lapply(probeset2entrezid,length)
temp.probeset<-names(probeset2entrezid)
temp.vect<-cbind(temp.probeset,as.character(temp.length))
temp.list<-apply(temp.vect,1,my.vect.construct)

temp1<-unlist(probeset2entrezid)
temp2<-unlist(temp.list)
entrezid2probeset<-cbind(temp1,temp2)
colnames(entrezid2probeset)<-c("entrezID","probeset")

####################
# generate the list that contains the probeset IDs for each pathway
pathway2probeset<-lapply(pathway2entrezid,my.list.construct,entrezid2probeset=entrezid2probeset)
names(pathway2probeset)<-sapply(names(pathway2probeset),my.pathway.namestand,replace="_")
msigdb.pathway.list<-pathway2probeset
return(msigdb.pathway.list)
}
#-----------------------------------------------------------------------------------------------------
# this function load in the probeset ids in each gene sets generated by chip2chip in GSEA
my.msgidb.list.load<-function(x){
temp<-readLines(x)
temp.names<-unname(sapply(temp,my.element.extract,splitchar="\t",index=1))
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp<-sapply(temp,my.element.remove,splitchar="\t",index=1)
temp.list<-sapply(temp,strsplit,split="\t")
names(temp.list)<-temp.names
return(temp.list)
}


#-----------------------------------------------------------------------------------------------------
# this function normalizes the gene expression data so that each gene has a 0 mean and 1 sd
my.normalize<-function(temp1){
temp<-apply(temp1,1,mean)
temp<-matrix(rep(temp,ncol(temp1)),nrow=nrow(temp1),byrow=F)
serum.data.norm<-temp1-temp
temp<-apply(serum.data.norm,1,sd)
temp<-matrix(rep(temp,ncol(temp1)),nrow=nrow(temp1),byrow=F)
serum.data.norm<-serum.data.norm/temp
return(serum.data.norm)
}


#-----------------------------------------------------------------------------------------------------
# this function applies the t-test on x by comparing x[1:n] to x[(n+1):length(x)]. The test will be two sided
my.tstat<-function(x,n){
temp<-t.test(x[1:n],x[(n+1):length(x)],alternative="two.sided")
return(tstat=temp$statistic)
}

my.ttest<-function(x,n){
if(sd(x)==0){
return(1)
}else{
temp<-t.test(x[1:n],x[(n+1):length(x)],alternative="two.sided")
return(temp$p.value)
}
}

#-----------------------------------------------------------------------------------------------------
# this function applies the t-test on x by comparing x[1:n] to x[(n+1):length(x)]. The test will be two sided
my.wilcox<-function(x,n){
temp<-wilcox.test(x[1:n],x[(n+1):length(x)],alternative="two.sided")
if(median(x[1:n])<median(x[(n+1):length(x)])){
result<--temp$p.value
}else{
result<-temp$p.value
}
return(result)
}


#-----------------------------------------------------------------------------------------------------
# this function applies the Kruskal-Wallis test on x by comparing groupings determined by clustering.result. 
# The test will be two sided
my.kruskal<-function(x,clustering.result){
return(kruskal.test(x,clustering.result,na.action="rm")$p.value)
}






#----------------------------------------------------------------------------------
# this function loads in the mclust results by the pathway based clustering method
my.result.load<-function(pathway.list,result.dir,compartment.prefix){
	###############
	# load in the results and calculate the classification matrix for sputum and blood
	sample.names<-character(0)
	mclust.list<-list()
	class.matrix<-numeric(0)
	pvalue.list<-list()
	minpvalue.vect<-numeric()
	pathway.names<-character(0)
	library(mclust)
	# load in the pvalues and extract the minimum pvalue for each valide pathway
	for(i in 1:length(pathway.list)){
		# load in the classification results
		filepath<-paste(result.dir,"/",compartment.prefix,"_",names(pathway.list)[i],"_mclust.txt",sep="")
		if(file.exists(filepath)==F){
			cat(filepath," does not exist!\n",sep="")
			next
		}
		temp<-readLines(filepath)
		if(temp=="NA"){	
			cat(filepath," is NA!\n",sep="")
			next	
		}
		temp<-matrix(unlist(strsplit(temp,split="\t")),nrow=length(temp),byrow=T)
		if(length(sample.names)==0){
			sample.names<-temp[,1]
		}else{
			if(length(setdiff(sample.names,temp[,1]))>0 | length(setdiff(temp[,1],sample.names))>0){
				cat(names(pathway.list)[i],": different sample names!\n")
				break
			}
		}
		rownames(temp)<-temp[,1]
		temp<-temp[sample.names,]
		mclust.list[[i]]<-temp
		class.matrix<-cbind(class.matrix,temp[,2])

		# load in the pvalues and extract the minimum p-value
#		filepath<-paste(result.dir,"/",compartment.prefix,"_",names(pathway.list)[i],".txt",sep="")
#		temp<-readLines(filepath)
#		temp<-matrix(unlist(strsplit(temp,split="\t")),nrow=length(temp),byrow=T)
#		temp.pvalue<-as.numeric(temp[,3])
#		temp.pvalue<-temp.pvalue[!is.na(temp.pvalue)]
#		minpvalue.vect<-c(minpvalue.vect,min(temp.pvalue))
#		pvalue.list[[length(pvalue.list)+1]]<-as.numeric(temp[,3])
		pvalue.list[[length(pvalue.list)+1]]<-NA
		pathway.names<-c(pathway.names,names(pathway.list)[i])
	}
	colnames(class.matrix)<-pathway.names
	names(pvalue.list)<-pathway.names
	
	return(list(class.matrix=class.matrix,pvalue.list=pvalue.list))
}






#----------------------------------------------------------------------------------
# this function calculates the distance matrix from the classmatrix generated by my.result.load
my.dist.matrix<-function(class.matrix,pathway.selected){
	# calculate the distances between samples using the classification results by selected pathways
	dist.matrix<-matrix(0,nrow=nrow(class.matrix),ncol=nrow(class.matrix))
	rownames(dist.matrix)<-rownames(class.matrix)
	colnames(dist.matrix)<-rownames(class.matrix)

	#pathway.selected<-pathway.names[substr(pathway.names,1,5)=="KEGG_"]
	for(i in 1:nrow(dist.matrix)){
		for(j in 1:nrow(dist.matrix)){
			sample1.name<-rownames(dist.matrix)[i]
			sample2.name<-rownames(dist.matrix)[j]
			temp1<-class.matrix[sample1.name,pathway.selected]
			temp2<-class.matrix[sample2.name,pathway.selected]
			dist.matrix[i,j]<-sum(temp1==temp2)
		}
	}
	return(dist.matrix)
}




#----------------------------------------------------------------------------------
# this function assess the significance of the correlation between a continous vector and a clustering results
# The two vectors have to be matched. 
# x can be a character vector but clustering result should be a numeric vector
my.clinic.con.pvalue<-function(x,clustering.result){
	library(caTools)
	cluster.num<-length(unique(clustering.result))
	cluster.vect<-unique(clustering.result)
#	result<-list(kw.pvalue=-2,anova.pvalue=-2,cluster.pvalue=rep(-2,max(as.numeric(clustering.result))),cluster.maxp=rep(-2,max(as.numeric(clustering.result))),cluster.minp=rep(-2,max(as.numeric(clustering.result))))
	result<-list(kw.pvalue=-2,anova.pvalue=-2,cluster.pvalue=rep(-2,cluster.num-1),cluster2cluster.pvalue=rep(-2,nrow(combs(1:cluster.num,2))),cluster2cluster.anovapvalue=rep(-2,nrow(combs(1:cluster.num,2))))

	temp<-as.numeric(x)
	if(sum(!is.na(temp))<5){
		return(result)
	}
	# calculate the p-value for the correlation based on kruskal-wallis test across all the given clusters (including controls if included)
#	temp.1<-temp[clustering.result!=0]
#	temp.2<-clustering.result[clustering.result!=0]
	temp.1<-temp
	temp.2<-clustering.result
	result$kw.pvalue<-kruskal.test(temp.1,temp.2,na.action="rm")$p.value
	# calculate the p-value for the correlation based on anova test
	if(0%in%temp.1){
	my.temp<-data.frame(temp1=log(temp.1+1e-06,base=2),temp2=as.factor(temp.2))
	}else{
	my.temp<-data.frame(temp1=log(temp.1,base=2),temp2=as.factor(temp.2))
	}

	temp.list<-split(my.temp[,1],my.temp[,2])
	temp.length<-unlist(lapply(temp.list,length))
	temp.list<-lapply(temp.list,is.na)
	temp.list<-unlist(lapply(temp.list,sum))
	temp.list<-temp.length-temp.list
	if(sum(temp.list<2)>0){
		result$anova.pvalue<-NA
	}else{
		my.fit<-oneway.test(temp1~temp2,data=my.temp)
		result$anova.pvalue<-my.fit$p.value
	}
	# calculate the pvalues by comparing each cluster with the combination of the other two clusters using wilcoxon test
#	cluster.num<-max(clustering.result)	
	for(j in 1:length(unique(temp.2))){
		# for each cluster, calculate the p-values assessing the difference between this cluster and all other clusters
		# we only focus on data with less than 10% missing data in every cluster
		if(sum(!is.na(temp.1[temp.2==unique(temp.2)[j]]))<5 | sum(!is.na(temp.1[temp.2!=unique(temp.2)[j]]))<5){
			result$cluster.pvalue[j]<-NA
			next
		}

		result$cluster.pvalue[j]<-wilcox.test(temp.1[temp.2==unique(temp.2)[j]],temp.1[temp.2!=unique(temp.2)[j]],na.action="rm")$p.value
		if(median(temp.1[temp.2==unique(temp.2)[j]],na.rm=T)<median(temp.1[temp.2!=unique(temp.2)[j]],na.rm=T)){
			result$cluster.pvalue[j]<--result$cluster.pvalue[j]
		}
	}
	names(result$cluster.pvalue)<-as.character(unique(temp.2))
	# calculate the 3 pvalues when comparing each cluster with each of the other two clusters
#	cluster.num<-max(clustering.result)
	library(caTools)
	temp.combs<-combs(1:cluster.num,2)
	temp.pvalue.vect<-numeric()
	temp.pvalue.direction<-numeric()
	for(j in 1:nrow(temp.combs)){
		# for each cluster, calculate the p-values assessing the difference between this cluster and each of the other two clusters
	
		# we only focus on data with less than 10% missing data in every cluster
		if(sum(!is.na(temp[clustering.result==cluster.vect[temp.combs[j,1]]]))<5 | sum(!is.na(temp[clustering.result==cluster.vect[temp.combs[j,2]]]))<5){
			temp.pvalue.vect<-c(temp.pvalue.vect,NA)
			temp.pvalue.direction<-c(temp.pvalue.direction,0)
			next
		}

		temp.pvalue.vect<-c(temp.pvalue.vect,wilcox.test(temp[clustering.result==cluster.vect[temp.combs[j,1]]],temp[clustering.result==cluster.vect[temp.combs[j,2]]],na.action="rm")$p.value)
		if(median(temp[clustering.result==cluster.vect[temp.combs[j,1]]],na.rm=T)<median(temp[clustering.result==cluster.vect[temp.combs[j,2]]],na.rm=T)){
			temp.pvalue.direction<-c(temp.pvalue.direction,-1)				
		}else{
			temp.pvalue.direction<-c(temp.pvalue.direction,1)				
		}
		#result$cluster.maxp[j]<-wilcox.test(temp[clustering.result==j],temp[clustering.result!=j],na.action="rm")$p.value
		#result$cluster.maxp[j]<-temp.pvalue.direction[temp.pvalue.vect==max(temp.pvalue.vect)]*max(temp.pvalue.vect)
		#result$cluster.minp[j]<-temp.pvalue.direction[temp.pvalue.vect==min(temp.pvalue.vect)]*min(temp.pvalue.vect)
		
	}
	result$cluster2cluster.pvalue<-temp.pvalue.vect*temp.pvalue.direction
	temp.names<-cbind(cluster.vect[temp.combs[,1]],cluster.vect[temp.combs[,2]])
	names(result$cluster2cluster.pvalue)<-apply(temp.names,1,paste,collapse="_")


	# use anova test to assess the difference between different clusters
	library(caTools)
	temp.combs<-combs(1:cluster.num,2)
	temp.pvalue.vect<-numeric()
	temp.pvalue.direction<-numeric()
	for(j in 1:nrow(temp.combs)){
		# for each cluster, calculate the p-values assessing the difference between this cluster and each of the other two clusters
	
		# we only focus on data with less than 10% missing data in every cluster
		if(sum(!is.na(temp[clustering.result==cluster.vect[temp.combs[j,1]]]))<5 | sum(!is.na(temp[clustering.result==cluster.vect[temp.combs[j,2]]]))<5){
			temp.pvalue.vect<-c(temp.pvalue.vect,NA)
			temp.pvalue.direction<-c(temp.pvalue.direction,0)
			next
		}

		temp.vect<-temp[clustering.result%in%c(cluster.vect[temp.combs[j,1]],cluster.vect[temp.combs[j,2]])]
		if(0%in%temp.vect){
		my.temp<-data.frame(temp1=log(temp.vect+1e-06,base=2),temp2=as.factor(clustering.result[clustering.result%in%c(cluster.vect[temp.combs[j,1]],cluster.vect[temp.combs[j,2]])]))
		}else{
		my.temp<-data.frame(temp1=log(temp.vect,base=2),temp2=as.factor(clustering.result[clustering.result%in%c(cluster.vect[temp.combs[j,1]],cluster.vect[temp.combs[j,2]])]))
		}

		temp.list<-split(my.temp[,1],my.temp[,2])
		temp.length<-unlist(lapply(temp.list,length))
		temp.list<-lapply(temp.list,is.na)
		temp.list<-unlist(lapply(temp.list,sum))
		temp.list<-temp.length-temp.list
		if(sum(temp.list<2)>0){
			temp.pvalue.vect<-c(temp.pvalue.vect,NA)
		}else{
			my.fit<-oneway.test(temp1~temp2,data=my.temp)
			temp.pvalue.vect<-c(temp.pvalue.vect,my.fit$p.value)
		}
		if(!is.na(temp.pvalue.vect[length(temp.pvalue.vect)])){
			if(median(temp[clustering.result==cluster.vect[temp.combs[j,1]]],na.rm=T)<median(temp[clustering.result==cluster.vect[temp.combs[j,2]]],na.rm=T)){
				temp.pvalue.direction<-c(temp.pvalue.direction,-1)				
			}else{
				temp.pvalue.direction<-c(temp.pvalue.direction,1)				
			}
		}
		#result$cluster.maxp[j]<-wilcox.test(temp[clustering.result==j],temp[clustering.result!=j],na.action="rm")$p.value
		#result$cluster.maxp[j]<-temp.pvalue.direction[temp.pvalue.vect==max(temp.pvalue.vect)]*max(temp.pvalue.vect)
		#result$cluster.minp[j]<-temp.pvalue.direction[temp.pvalue.vect==min(temp.pvalue.vect)]*min(temp.pvalue.vect)
		
	}
	result$cluster2cluster.anovapvalue<-temp.pvalue.vect*temp.pvalue.direction
	temp.names<-cbind(cluster.vect[temp.combs[,1]],cluster.vect[temp.combs[,2]])
	names(result$cluster2cluster.anovapvalue)<-apply(temp.names,1,paste,collapse="_")

	return(result)
}




#----------------------------------------------------------------------------------
# this function assess the significance of the correlation between a categorical vector and a clustering results
# The two vectors have to be matched. 
# x can be a character vector but clustering result should be a numeric vector
my.clinic.cat.pvalue<-function(x,clustering.result){
	result<-list(pvalue=-2,trend.test.pvalue=-2,cluster.pvalue=rep(-2,max(as.numeric(clustering.result))))

	temp<-as.numeric(as.factor(x))
	if(sum(!is.na(temp))<5){
		return(result)
	}
	# calculate the p-value for the correlation based on chi square test
	result$pvalue<-chisq.test(temp,clustering.result)$p.value
	if(length(as.numeric(table(x)))==2){
		temp1<-as.matrix(table(x,clustering.result))
		temp1<-temp1[,order(temp1[1,]/apply(temp1,2,sum),decreasing=T)]
		result$trend.test.pvalue<-prop.trend.test(temp1[1,],apply(temp1,2,sum))$p.value
	}
	return(result)
}


#----------------------------------------------------------------------------------
# this function extract the element at position index in a given list
my.list.element.extract<-function(x,index=1){
if(length(x)==0){
cat("Given list is empty!\n")
return(NA)
}
if(index>length(x)){
cat("Not enough elements in X!\n")
return(NA)
}

if(index<0){
return(x[[length(x)]])
}else{
return(x[[index]])
}


}

#----------------------------------------------------------------------------------
# This function generate colors for categorical variables
my.color.gen.cat<-function(date){
	all.color.vect = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
	all.color.vect = all.color.vect[grep('white', all.color.vect, invert = T)]
	
	if(length(unique(date))>length(all.color.vect)){
			cat(length(all.color.vect)," does not have enough colors for the variable!\n",sep="")
			return(NA)
	}
	
	library(gplots)
	# this function generates the color names for every element in date
	temp<-sort(unique(date),decreasing=T)
	color.pool<-sample(all.color.vect,length(temp))
	color.int<-rep(0,length(date))
	color.num<-rep("",length(date))
	
	for(i in 1:length(temp)){
		#color.num[date==temp[i]]<-colorpanel(length(temp),low="red",mid="grey",high="blue")[i]
		color.num[date==temp[i]]<-color.pool[i]
		color.int[date==temp[i]]<-i
	}
	color.names<-color.pool
	result<-list(color.num=color.num,color.int=color.int,color.names=color.names)
	return(result)
}

#----------------------------------------------------------------------------------
# this function generate a continuous spectrum of colors from blue to red for values in the variable
my.color.gen.con<-function(date,nbreaks=500,method="absolute"){
# method represent how the colors would be assigned. absolute means it represent the real number while quantile means it represent the quantile instead of the original values
library(gplots)

if(method=="absolute"){
	# generate the intervals and the corresponding colors
	color.spectrum<-colorpanel(nbreaks,low="blue",mid="grey",high="red")
	interval.length<-(max(date,na.rm=T)-min(date,na.rm=T))/(nbreaks-1)
	
	# generate the color for each element in date
	interval.index<-sapply(date,function(x,ilength,imin){return(ceiling((x-imin+0.5*ilength)/ilength))},ilength=interval.length,imin=min(date,na.rm=T))
	interval.index[!is.na(interval.index) & interval.index==0]<-1
	color.num<-color.spectrum[interval.index]
	color.names<-color.num
	color.int<-interval.index
	result<-list(color.num=color.num,color.int=color.int,color.names=color.names)
	return(result)
}else{
		
	# this function generates the color names for every element in date
	temp<-sort(date,decreasing=F)
	color.int<-rep(0,length(date))
	color.num<-rep("",length(date))
	for(i in 1:length(temp)){
		color.num[date==temp[i]]<-colorpanel(length(temp),low="blue",mid="grey",high="red")[i]
		color.int[date==temp[i]]<-i
	}
	color.names<-temp
	result<-list(color.num=color.num,color.int=color.int,color.names=color.names)
	return(result)
		
}
	
}

















