source("temp0.0.R")

genes=read.table("genes.txt",sep="\t")
allgene.matrix=matrix(0,ncol=3)
allgene.matrix=allgene.matrix[-1,]
stime=Sys.time()
for (k in 1:(dim(genes))[1])
 {
full_matrix=matrix(0,ncol=50)
full_matrix=full_matrix[-1,]
colnames(full_matrix)=c(paste(c("Name","SNP","Chr","MapInfo","cpg.name"),sep=""),1:45)


     #k=as.numeric(b)
     cat("gene.no =",k,"\n")
     # gene name
     #Gene.name="TSC1"
     Gene.name=as.character(genes[k,])
     
     #cat("Hello8_",Gene.name,"\n")
     
     # gene expr
     exprn=matrix(NA,nrow=length(mat.id),ncol=1)
     exprn[match(colnames(case.norm.count.50perc.miss)[-(1:4)],mat.id)]=round(as.numeric(case.norm.count.50perc.miss[which(case.norm.count.50perc.miss[,1]==Gene.name)[1],-(1:4)]),4)
     exprn[match(colnames(control.norm.count.nomiss)[-(1:4)],mat.id)]=round(as.numeric(control.norm.count.nomiss[which(control.norm.count.nomiss[,1]==Gene.name)[1],-(1:4)]),4)
     
     # standardize the exprn values (based on non missing gene expression) accross case and control individuals
     # (exprn-min(exprn))/max(exprn)  is case of no missing
     
     avail.exprn=exprn[which(is.na(exprn)==FALSE)]
     #avail.exprn=log(avail.exprn)
     exprn1=matrix(NA,nrow=dim(exprn)[1],ncol=1)
     exprn1[which(is.na(exprn)==FALSE)]=(avail.exprn-mean(avail.exprn))/sd(avail.exprn)
     exprn=exprn1

exprn=c(rep(-99,5),exprn)
full_matrix=rbind(full_matrix,t(exprn))


     # Find the snps in 2000 bp upstream and 2000 bp downstream of the gene start and end. If
     store.chr=NULL
     store.chr=as.character(control.norm.count.nomiss[which(control.norm.count.nomiss[,1]==Gene.name),2])
     store.gene.start=control.norm.count.nomiss[which(control.norm.count.nomiss[,1]==Gene.name),3]-2000
     store.gene.end=control.norm.count.nomiss[which(control.norm.count.nomiss[,1]==Gene.name),4]+2000
     
     
     case.geno.file.chr=case.geno.file[which(case.geno.file[,3]==(unlist(strsplit(store.chr,split="chr"))[2])),]
     control.geno.file.chr=control.geno.file[which(control.geno.file[,3]==(unlist(strsplit(store.chr,split="chr"))[2])),]
     
     id.start=which(control.geno.file.chr[,4]>store.gene.start[1])
     id.end=which(control.geno.file.chr[,4]<store.gene.end[1])
     
     id=as.numeric(na.omit(match(id.start,id.end)))#which(duplicated(c(id.start,id.end))==TRUE)
     
     case.geno.file.gene=case.geno.file.chr[id,]
     control.geno.file.gene=control.geno.file.chr[id,]

#cc=cbind(control.geno.file.gene,case.geno.file.gene[,-(1:4)])
#full_matrix=rbind(full_matrix,cc)
 
     
     
     # We start with gene TSC1 in Chr 9
     
     func_meth=function(jj)
     {
         methyl=matrix(NA,nrow=length(mat.id),ncol=1)
         methyl[match(colnames(meth.control)[-c(1,9:11)],mat.id)]=(meth.control[which(meth.control[,9]==Gene.name)[jj],-c(1,9:11)])
         methyl[match(colnames(meth.case)[-c(1,23:25)],mat.id)]=(meth.case[which(meth.case[,23]==Gene.name)[jj],-c(1,23:25)])
         return(methyl)
         
     }
     # genotype
     
     
     func_geno=function(ii)
     {
         geno=matrix(NA,nrow=length(mat.id),ncol=1)
         geno[as.numeric(na.omit(match(colnames(case.geno.file.gene),mat.id)))]=as.numeric(case.geno.file.gene[which(case.geno.file.gene[,3]==(unlist(strsplit(store.chr,split="chr"))[2]))[ii],-(1:4)])
         geno[as.numeric(na.omit(match(colnames(control.geno.file.gene),mat.id)))]=as.numeric(control.geno.file.gene[which(control.geno.file.gene[,3]==(unlist(strsplit(store.chr,split="chr"))[2]))[ii],-(1:4)])
         return(geno)
         
     }
     x.cpgs=methyl_all.temp=NULL
     x.cpgs=1:length(which(meth.control[,9]==Gene.name))
     methyl_all.temp=do.call("cbind",lapply(x.cpgs,func_meth))
     
     func.methval.std=function(x)
     {
         methval=unlist(methyl_all.temp[,x])
         avail.meth=methval[which(is.na(methval)==FALSE)]
         avail.meth1=matrix(NA,nrow=length(methval),ncol=1)
         avail.meth1[which(is.na(methval)==FALSE)]=(avail.meth-mean(avail.meth))/sd(avail.meth)
         return(avail.meth1)
     }
     
     y.geno=cpgs=methyl_all=geno_all=NULL
     temp1=temp2=NULL
     cpg.name=geno_all_temp2=store_methnames=methyl_all_temp1=NULL
     cpg.name=rep(-99,dim(case.geno.file.gene)[1])


     
     
     if (dim(case.geno.file.gene)[1]>=1)
     {
         y.geno=1:(dim(case.geno.file.gene)[1])
         geno_all=do.call("cbind",lapply(y.geno,func_geno))
         temp2=1:length(y.geno)
         geno_all_temp2=cbind(case.geno.file.gene[,1:4],cpg.name,t(geno_all))
         colnames(geno_all_temp2)=c(paste(c("Name","SNP","Chr","MapInfo","cpg.name"),sep=""),1:45)
         full_matrix=rbind(full_matrix,geno_all_temp2)

     }
     if (dim(case.geno.file.gene)[1]<1)  y.geno=-99
     if (dim(methyl_all.temp)[2]>=1)
     {
         cpgs=1:(dim(methyl_all.temp)[2])
         methyl_all=do.call("cbind",lapply(cpgs,func.methval.std))
         temp1=1:length(cpgs)
         store_methnames=meth.control[which(meth.control[,9]==Gene.name),1]
         methyl_all_temp1=data.frame(matrix(-99,ncol=4,nrow=length(store_methnames)),store_methnames,t(methyl_all))
         colnames(methyl_all_temp1)=c(paste(c("Name","SNP","Chr","MapInfo","cpg.name"),sep=""),1:45)
         full_matrix=rbind(full_matrix,methyl_all_temp1)

     }
     if (dim(methyl_all.temp)[2]<1)  cpgs=-99

fname=paste(c("geneid_",k,"_",Gene.name,".txt"),collapse="")
write.table(full_matrix,fname,quote=F,sep="\t",row.names=F,col.names=F)


if (length(y.geno)>1 && length(cpgs)>1)
allgene.matrix=rbind(allgene.matrix,c(k,length(y.geno),length(cpgs)))

 
if (length(y.geno)==1 && y.geno==-99 && length(cpgs)>=1 && cpgs!=-99)
allgene.matrix=rbind(allgene.matrix,c(k,0,length(cpgs)))

if (length(y.geno)==1 && y.geno!=-99 && length(cpgs)>=1 && cpgs!=-99)
allgene.matrix=rbind(allgene.matrix,c(k,1,length(cpgs)))

if (length(y.geno)>=1 && y.geno!=-99 && length(cpgs)==1 && cpgs!=-99)
allgene.matrix=rbind(allgene.matrix,c(k,length(y.geno),1))

if (length(y.geno)>=1 && y.geno!=-99 && length(cpgs)==1 && cpgs==-99)
allgene.matrix=rbind(allgene.matrix,c(k,length(y.geno),0))

if (length(y.geno)==1 && y.geno==-99 && length(cpgs)==1 && cpgs==-99)
allgene.matrix=rbind(allgene.matrix,c(k,0,0))


write.table(allgene.matrix,"allgene.matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)
 }



genes1=genes[allgene.matrix[,1],]
allgene.matrix=as.data.frame(allgene.matrix,genes1)
allgene.matrix=cbind(allgene.matrix[,1],genes1,allgene.matrix[,-1])
colnames(allgene.matrix)=paste(c("GeneID","Gene.name","num.Genotype","num.Methylation"))

write.table(allgene.matrix,"allgene.matrix.txt",quote=F,sep="\t",row.names=F,col.names=T)

Sys.time()-stime

