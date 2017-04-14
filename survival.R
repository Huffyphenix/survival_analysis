#!/usr/bin/Rscript
############################ get exp data ############################
args<-commandArgs(TRUE)
exp<-read.table(args[1],sep="\t",row.names = 1,header = T)
disease_samples<-grep("0.",substr(colnames(exp),14,15))
disease_exp<-exp[,disease_samples]
colnames(disease_exp)<-substr(colnames(disease_exp),9,12)
disease_exp<-rbind(colnames(disease_exp),disease_exp)
rownames(disease_exp)[1]<-c("sampleID")
rem <- function(x,y){
  x <- as.matrix(x)
  r <- as.numeric(apply(x,y,function(i) sum(i == 0)))
  if(y==1){remove <- which(r > dim(x)[2]*0.5)}
  else{remove <- which(r == dim(x)[1])}
  return(remove)
}
del_1<-rem(disease_exp,1)
if(length(del_1)!=0){disease_exp<-disease_exp[-del_1,]}
del_2<-rem(disease_exp,2)
if(length(del_2)!=0){disease_exp<-disease_exp[,-del_1]}
############################ get clinical data ############################
survival<-read.table(args[2],row.names = 1,quote="",na.strings = "NA",sep="\t")
del<-as.numeric(which(apply(survival,2,function(x)all(is.na(x)))=="TRUE"))
if(length(del)!=0){survival<-survival[,-del]}
del<-as.numeric(which(is.na(survival["patient.vital_status",])=="TRUE"))
if(length(del)!=0){survival<-survival[,-del]}
patient_id<-as.vector(apply(survival["patient.patient_id",],2,function(i) toupper(i)))
patient.days_to_death<-matrix(apply(survival["patient.days_to_death",],2,function(i) as.numeric(i)),nrow=1)
patient.days_to_last_followup<-matrix(apply(survival["patient.days_to_last_followup",],2,function(i) as.numeric(i)),nrow=1)
surtime<-function(death,follow_up){
  time<-matrix(NA,nrow = 1,ncol=ncol(death))
  for(i in 1:ncol(death)){
    if(!is.na(death[1,i])){time[1,i]=death[1,i]}
    else{time[1,i]=follow_up[1,i]}
  }
  return(time)
}
time<-as.numeric(round(surtime(patient.days_to_death,patient.days_to_last_followup)/30,3))
status<-as.vector(apply(survival["patient.vital_status",],2,function(i){if(i=="alive") return(0) else return(1)}))
survival_clinical<-data.frame(sampleID=patient_id,OS=as.numeric(time),status=status)
del<-as.numeric(which(is.na(survival_clinical[,"OS"])=="TRUE"))
if(length(del)!=0){survival_clinical<-survival_clinical[-del,]}
del<-as.numeric(which(survival_clinical[,"OS"]<=0))
if(length(del)!=0){survival_clinical<-survival_clinical[-del,]}
############################ merge exp and clinical ############################
OS_exp<-merge(survival_clinical,t(disease_exp),by="sampleID")
del=c()
for(j in 4:ncol(OS_exp)){
  le=length(as.numeric(which(is.na(OS_exp[,j])=="TRUE")))
  D=nrow(OS_exp)-le
  if(D<=5){del<-c(del,j)}
}
if(length(del)!=0){OS_exp<-OS_exp[,-del]}

del=c()
for(j in 1:nrow(OS_exp)){
  le=length(as.numeric(which(is.na(OS_exp[j,])=="TRUE")))
  D=ncol(OS_exp)-le
  if(D==0){del<-c(del,j)}
}
if(length(del)!=0){OS_exp<-OS_exp[,-del]}

############################ do survival analysis ############################
p_signif<-function(p){
  if(p>0.05) label=c("")
  if(p>0.05&&p<=0.1) label=c(".")
  if(p<=0.05&&p>0.01) label=c("*")
  if(p<=0.01&&p>0.001) label=c("**")
  if(p<=0.001&&p>0.0001) label=c("***")
  if(p<=0.0001&&p>0.00001) label=c("****")
  if(p<=0.00001) label=c("*****")
  return(label)
}
library("survival")
Result<-matrix(NA,nrow=ncol(OS_exp)-3,ncol=13,byrow = TRUE)
rownames(Result) <- colnames(OS_exp)[4:ncol(OS_exp)]
colnames(Result) <- c("KMp", "N","low_N","high_N","low_exp","high_exp","median_exp","log2FC","prognosis","low_sur_median","high_sur_median","low_surv.rate","high_surv.rate")
p_signif_list<-c()
boundary<-as.numeric(as.character(args[6]))
pvalue<-as.numeric(as.character(args[5]))
for(j in 4:ncol(OS_exp)){
  tmp_exp<-as.numeric(as.character(OS_exp[,j]))
  tmp_OS<-OS_exp$OS
  tmp_status<-OS_exp$status
  tmp_data<-data.frame(OS=tmp_OS,status=tmp_status,exp=tmp_exp)
  del<-which(is.na(tmp_data$exp))
  if(length(del)!=0 & length(del)<length(tmp_exp)){tmp_data[-del,]->tmp_data}
  tmp_group<-as.factor(ifelse(as.numeric(as.character(tmp_data$exp)) <= quantile(as.numeric(as.character(tmp_data$exp)),boundary), "Low","High"))
  tmp_data<-data.frame(tmp_data,group=as.factor(tmp_group))
  surv<-Surv(tmp_data$OS,tmp_data$status)
  #tempresult<-try(c<-coxph(surv ~ exp,tmp_data),silent=TRUE)
  i=j-3
  #if(!inherits(tempresult, "try-error")){
#  Result[i, c("coef", "Exp(coef)", "Coxp")] <- summary(c)$coefficients[1,c("coef", "exp(coef)", "Pr(>|z|)" )]
  model1 <- survdiff(surv~group,data=tmp_data,na.action = na.exclude)
  Result[i, c("KMp")] <- signif(1-pchisq(model1$chisq, df=length(levels(factor(tmp_group)))-1),3)
  Result[i, "low_exp"] <- round(as.numeric(as.character(mean(tmp_data[which(tmp_data$group == "Low"),"exp"]))),2)
  Result[i, "high_exp"] <- round(as.numeric(as.character(mean(tmp_data[which(tmp_data$group == "High"),"exp"]))),2)
  Result[i, "log2FC"] <- round(log2((as.numeric(as.character(Result[i,"high_exp"]))+1)/(as.numeric(as.character(Result[i,"low_exp"]))+1)),3)
  fit <- survfit(surv~group, data=tmp_data, na.action=na.exclude)
  Result[i, "low_N"]<-fit$n[1]
  Result[i, "high_N"]<-fit$n[2]
  Result[i, "N"] <- sum(fit$n)
  Result[i,"median_exp"] <- round(median(tmp_data[,"exp"]),3)
  s.fit<-summary(fit)
  fit.data<-cbind(median=s.fit$table[,"median"],events=s.fit$table[,"events"],n=s.fit$table[,"n.start"])
  survrate<-(fit.data[,"n"]-fit.data[,"events"])/fit.data[,"n"]
  fit.data<-cbind(fit.data,surv.rate=survrate)
  if(is.na(fit.data["group=High","median"]) | is.na(fit.data["group=Low","median"])){
    if(fit.data["group=High","surv.rate"]<fit.data["group=Low","surv.rate"] ){Result[i, "prognosis"]="poor"}
    if(fit.data["group=High","surv.rate"]>fit.data["group=Low","surv.rate"] ){Result[i, "prognosis"]="good"}}else{
    if(fit.data["group=High","median"]<fit.data["group=Low","median"] & fit.data["group=High","surv.rate"]<fit.data["group=Low","surv.rate"] ){Result[i, "prognosis"]="poor"}
    if(fit.data["group=High","median"]>fit.data["group=Low","median"] & fit.data["group=High","surv.rate"]>fit.data["group=Low","surv.rate"] ){Result[i, "prognosis"]="good"}
#    if(fit.data["group=High","median"]<fit.data["group=Low","median"]){Result[i, "prognosis"]="poor"}
#    if(fit.data["group=High","median"]>fit.data["group=Low","median"]){Result[i, "prognosis"]="good"}
  }
  Result[i,"low_sur_median"]=fit.data["group=Low","median"]
  Result[i,"high_sur_median"]=fit.data["group=High","median"]
  Result[i,"low_surv.rate"]=round(fit.data["group=Low","surv.rate"],3)
  Result[i,"high_surv.rate"]=round(fit.data["group=High","surv.rate"],3)
  if(as.numeric(as.character(Result[i, c("KMp")])) <= pvalue){
    if(!is.na(Result[i,"prognosis"])){
    if(Result[i, "prognosis"]=="good"){
    pdf_name<-paste(args[7],"/",args[3],"/","PDF/good_prognosis","/",args[3],"_",args[4],"_",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][1],".",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][2],".pdf",sep="")}else{
    pdf_name<-paste(args[7],"/",args[3],"/","PDF/poor_prognosis","/",args[3],"_",args[4],"_",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][1],".",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][2],".pdf",sep="")}}else{
    pdf_name<-paste(args[7],"/",args[3],"/","PDF","/",args[3],"_",args[4],"_",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][1],".",strsplit(colnames(OS_exp)[j],split = "\\|")[[1]][2],".pdf",sep="")}
    pdf(pdf_name,width=8,height = 6)
    p_signif_list<-c(p_signif_list,colnames(OS_exp)[j])
    plot(fit,col=c("red","green"),lty=1,lwd=2,mark.time=TRUE,
         main=paste("Kaplan-Meier Curves of ",rownames(Result)[i]," in ",args[3],sep=""),
         xlab = "Survival in months",ylab="Survival rate",cex.lab=1.3,cex.axis=1.2)
    legend("bottomleft",c(paste(levels(as.factor(tmp_data$group))," n=",fit$n,sep="")),
                        col=c("red","green"),lty=c(1,1,0),bty="n",lwd=3)
    text(max(tmp_data$OS)-25,1,paste("logrank P =",signif(as.numeric(Result[i, c("KMp")]),digits = 2),p_signif(as.numeric(Result[i, c("KMp")])),sep=" "),cex=1)
    dev.off()
}}
##################write pvalue ordered RESULT##########################
Result<-Result[order(as.numeric(as.character(Result[,1]))),]
result.name<-paste(args[7],"/",args[3],"/",args[3],"_",args[4],"_survival.info.txt",sep="")
write.table(Result,result.name,quote=F,sep="\t")
#write pvalue ordered p_signif gene info
p_signif_gene_info<-Result[p_signif_list,]
p_signif_gene_info<-p_signif_gene_info[order(as.numeric(as.character(p_signif_gene_info[,1]))),]
psignif_list_name<-paste(args[7],"/",args[3],"/",args[3],"_",args[4],"_survival.p_signif.info.txt",sep="")
write.table(p_signif_gene_info,psignif_list_name,quote=F,sep="\t")
####################write pvalue ordered p_signif gene list#################
p_signif_gene_list<-rownames(p_signif_gene_info)
psignif_list_name<-paste(args[7],"/",args[3],"/",args[3],"_",args[4],"_survival.p_signif.list.txt",sep="")
write.table(p_signif_gene_list,psignif_list_name,quote=F,sep="\t",row.names=F)
#exp filter and FC filter
FC=1
p_signif_gene_info[as.numeric(as.character(p_signif_gene_info[,"log2FC"]))>=1,]->p_signif_gene_info
if(length(grep("mir",args[4],ignore.case=TRUE))!=0){
   exp=50
   p_signif_gene_info[as.numeric(as.character(p_signif_gene_info[,"high_exp"]))>=exp,]->p_signif_gene_info
}else if(length(grep("^rna$",args[4],ignore.case=TRUE))!=0){
   exp=1000
   p_signif_gene_info[as.numeric(as.character(p_signif_gene_info[,"high_exp"]))>=exp,]->p_signif_gene_info
}
psignif_list_name<-paste(args[7],"/",args[3],"/",args[3],"_",args[4],"_survival.p_signif","_exp",exp,"_FC",2*FC,".info.txt",sep="")
write.table(p_signif_gene_info,psignif_list_name,quote=F,sep="\t")
