#!/usr/bin/Rscript
args<-commandArgs(TRUE)

mirna2gene<-try(read.table(args[1],sep="\t"),silent=TRUE)
TF2gene<-try(read.table(args[2],sep="\t"),silent=TRUE)
miRNA2TF<-try(read.table(args[3],sep="\t"),silent=TRUE)
TF2mirna<-try(read.table(args[4],sep="\t"),silent=TRUE)
if(inherits(mirna2gene, "try-error")){mirna2gene<-matrix(ncol=3)}
if(inherits(TF2gene, "try-error")){TF2gene<-matrix(ncol=3)}
if(inherits(miRNA2TF, "try-error")){miRNA2TF<-matrix(ncol=3)}
if(inherits(TF2mirna, "try-error")){TF2mirna<-matrix(ncol=3)}
p_signif_info<-try(read.table(args[5],sep="\t",row.names=1,header = T),silent=TRUE)
if(inherits(p_signif_info, "try-error")){p_signif_info<-matrix(ncol=3)}
target_info<-matrix(NA,nrow = nrow(p_signif_info),ncol=4,byrow=TRUE)
rownames(target_info)=rownames(p_signif_info)
colnames(target_info)<-c("targeted_by_mirna","targeted_by_TF","targets","gene_id")
rbind(mirna2gene,miRNA2TF)->mirna_targets
rbind(TF2gene,TF2mirna)->TF_targets
rbind(mirna_targets,TF_targets)->gene_targets
#RNA targets
for(j in rownames(p_signif_info)){
  i=strsplit(j,'[|]')[[1]][1]
  target_info[j,"gene_id"]=i
  if(paste0(mirna_targets[grep(TRUE,mirna_targets[,2]==i),1],collapse = ",")!=""){
  target_info[j,"targeted_by_mirna"]<-paste0(mirna_targets[grep(TRUE,mirna_targets[,2]==i),1],collapse = ",")}
  if(paste0(TF_targets[grep(TRUE,TF_targets[,2]==i),1],collapse = ",")!=""){
  target_info[j,"targeted_by_TF"]<-paste0(TF_targets[grep(TRUE,TF_targets[,2]==i),1],collapse = ",")}
  if(paste0(gene_targets[grep(TRUE,gene_targets[,1]==i),2],collapse = ",")!=""){
  target_info[j,"targets"]<-paste0(gene_targets[grep(TRUE,gene_targets[,1]==i),2],collapse = ",")}
}
cbind(gene_id=target_info[,"gene_id"],p_signif_info)->p_signif_info
cbind(p_signif_info,target_info[,1:3])->p_signif_info
psignif_name<-paste(args[6],"/",args[7],"/",args[7],"_RNA_survival.p_signif.target_info.txt",sep="")
write.table(p_signif_info,psignif_name,quote=F,sep="\t")
#miRNA targets
p_signif_info<-try(read.table(args[8],sep="\t",row.names=1,header = T),silent=TRUE)
if(inherits(p_signif_info, "try-error")){p_signif_info<-matrix(ncol=3)}
target_info<-matrix(NA,nrow = nrow(p_signif_info),ncol=3,byrow=TRUE)
rownames(target_info)=rownames(p_signif_info)
colnames(target_info)<-c("targeted_by_TF","targets","gene_id")
for(j in rownames(p_signif_info)){
  i=strsplit(j,'[|]')[[1]][1]
  target_info[j,"gene_id"]=i
  if(paste0(mirna_targets[grep(TRUE,mirna_targets[,1]==i),2],collapse = ",")!=""){
  target_info[j,"targets"]<-paste0(mirna_targets[grep(TRUE,mirna_targets[,1]==i),2],collapse = ",")}
  if(paste0(TF_targets[grep(TRUE,TF_targets[,2]==i),1],collapse = ",")!=""){
  target_info[j,"targeted_by_TF"]<-paste0(TF_targets[grep(TRUE,TF_targets[,2]==i),1],collapse = ",")}
}
cbind(gene_id=target_info[,"gene_id"],p_signif_info)->p_signif_info
cbind(p_signif_info,target_info[,1:2])->p_signif_info
psignif_name<-paste(args[6],"/",args[7],"/",args[7],"_miRNA_survival.p_signif.target_info.txt",sep="")
write.table(p_signif_info,psignif_name,quote=F,sep="\t")



