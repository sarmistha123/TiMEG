# The following code is for TSC Brain samples only.

control.geno.file=read.table("control.geno.file.txt",sep="\t",header=T)
case.geno.file=read.table("case.geno.file.txt",sep="\t",header=T)

# Gene expressions for case
case.norm.count.50perc.miss=read.table("case.log.count.50perc.miss.txt",sep="\t",header=T)
J=dim(case.norm.count.50perc.miss)[2]-4
J1=unlist(strsplit(colnames(case.norm.count.50perc.miss)[-(1:4)],split="Subject_ID_"))[seq(2,J*2,by=2)]
o=order(as.numeric(J1))
case.norm.count.50perc.miss=case.norm.count.50perc.miss[,c(1:4,o+4)]
# We remove Subjects 34, 82 as there is no genotype information for these Subjects
case.norm.count.50perc.miss=case.norm.count.50perc.miss[,-c(13,25)]

# Gene expressions for control
control.norm.count.nomiss=read.table("control.log.count.nomiss.txt",sep="\t",header=T)
control.norm.count.nomiss=control.norm.count.nomiss[,c(1:4,10,8,9,7,6,11,5)]

# Methylations for case

meth.case=read.table("meth.case.txt",sep="\t",header=T)
# We remove Subject 34 as there is no genotype information for this Subject
meth.case=meth.case[,-12]
# CpGsites are till row 478524 in the methylation file
meth.case=meth.case[1:dim(meth.case)[1],]
meth.case=cbind(rownames(meth.case),meth.case)
colnames(meth.case)[1]=paste("CpG")

# Methylations for control

meth.control=read.table("meth.control.txt",sep="\t",header=T)
# CpGsites are till row 478524 in the methylation file
meth.control=meth.control[1:dim(meth.control)[1],]

meth.control=cbind(rownames(meth.control),meth.control)
colnames(meth.control)[1]=paste("CpG")

# Read phenotype, age and gender
covariates_file=read.table("covariates_file.txt",sep="\t",header=T)
# replace 0 disease status with -1
covariates_file[which(covariates_file[,6]==0),6]=-1

mat.id=NULL
mat.id=c(colnames(control.geno.file)[-(1:4)],colnames(case.geno.file)[-(1:4)])
write.table(mat.id,"mat.id.txt",sep="\t",quote=F,row.names=F,col.names=F)

# covariates and phenotype/disease status (Control: -1, Case: 1)
cov=covariates_file[match(mat.id,paste("Subject_ID_",covariates_file[,2],sep="")),c(6,5,4)]

#Gene.name="TSC1"

genes=meth.control[na.omit(match(control.norm.count.nomiss[,1],meth.control[,9])),9]

write.table(genes,"genes.txt",quote=F,sep="\t",col.names=F,row.names=F)

