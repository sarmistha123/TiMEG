args = commandArgs(trailingOnly=TRUE)
b<-args
print(b)
bb<-read.table(paste(as.character(b)),header=F)
#print(dim(bb))
fout<-paste("result_",b,sep="")
##### Write the main code below
### Start main code


genesymbol=unlist(strsplit(b,split=".",fixed=T))[1]
#print(genesymbol)
############################  Likelihood analysis code starts  ####################
f=function(X,data_matrix)
{
    bzg=X[1];
    bg=X[2];bm=X[3];be=X[4];
    ag=X[5];
    gg=X[6];gm=X[7];
    s1_2=X[8];
    
    s2_2=X[9];
    
    p=X[10];a0=X[11];g0=X[12];b0=X[13]
    if(p==0) p=p+0.00001
    if(p==1) p=p-0.00001
    
    if(s1_2==0) s1_2=s1_2+0.00001
    if(s2_2==0) s2_2=s2_2+0.00001
    
    A=NULL
    
    G=NULL; A1=NULL
    
    G=data_matrix[,3]
    
    A1= G*log(p) + (2-G)*log(1-p)
    A=A1
    
    # only M missing
    Y=E=G=M=Z_gender=NULL;A2=NULL
    a2=NULL
    a2=which(is.na(data_matrix[,4]=="NA")==TRUE)[which(is.na( match( which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]))==TRUE)]
    
    
    if(length(a2)==0)    A=A
    if(length(a2)>0)
    {
        Y=data_matrix[a2,1]
        G=data_matrix[a2,3]
        E=data_matrix[a2,5]
        Z_gender=data_matrix[a2,2]
        
        A2=-log(1+exp(-(((Y*pi/sqrt(3))*(E*(bm * gm*s1_2 + be*s2_2 + be*s1_2*gm^2)+(s2_2 + gm^2*s1_2)*(b0 +  bzg*Z_gender + bg*G)+ bm*s2_2*(a0 + ag*G)- bm*gm*s1_2*(g0 + gg*G)))/(sqrt(s2_2 + gm^2*s1_2)*sqrt(bm^2*s1_2*s2_2 + (pi*pi/3)*(s2_2 + gm^2*s1_2)))))) - (((1/2)*(E-(g0 + gg*G + gm*(a0 + ag*G)))^2)/(s2_2+(s1_2)*(gm^2)))-(1/2)*log(2*pi*(s2_2+(gm^2)*(s1_2)))
        
        A[a2] = A[a2]+A2
        
    }
    
    
    # only E missing
    Y=G=M=Z_gender=NULL;A3=NULL
    a3=NULL
    a3=which(is.na(data_matrix[,5]=="NA")==TRUE)[which(is.na( match( which(is.na(data_matrix[,5]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]))==TRUE)]
    
    if(length(a3)==0)    A=A
    if(length(a3)>0)
    {
        
        Y=data_matrix[a3,1]
        G=data_matrix[a3,3]
        M=data_matrix[a3,4]
        Z_gender=data_matrix[a3,2]
        
        A3=-log(1+exp(-((Y*(pi/sqrt(3))*(b0  +  bzg*Z_gender + bg*G + bm*M + be*(g0 + gg*G + gm*M)))/sqrt((pi*pi/3) + (be^2) * s2_2)))) - (1/2)*(((M - a0 - ag*G)^2)/s1_2) - (1/2)*log(2*pi*s1_2)
        
        A[a3] = A[a3]+A3
    }
    
    
    
    # none of E or M missing
    Y=E=G=M=Z_gender=NULL;A4=NULL
    a4=NULL
    a4=-sort(unique(c(which(is.na(data_matrix[,5])=="TRUE"),which(is.na(data_matrix[,4])=="TRUE"))))
    
    if(length(which(is.na(data_matrix[,5])=="TRUE"))==0 && length(which(is.na(data_matrix[,4])=="TRUE"))==0)
    {
        for(i in 1:length(data_matrix[,5]))
        a4=c(a4,i)
    }
    #length(data_matrix[,6])
    
    
    if(length(a4)==0)    A=A
    if(length(a4)>0)
    {
        Y=data_matrix[a4,1]
        G=data_matrix[a4,3]
        M=data_matrix[a4,4]
        E=data_matrix[a4,5]
        Z_gender=data_matrix[a4,2]
        
        A4=  -log(1+exp(-(Y*(b0  + bzg*Z_gender + bg*G + bm*M + be*E)))) - (1/2)*(((E- g0 - gm*G - gm*M)^2)/s2_2) -(1/2)*log(2*pi*s2_2) - (1/2)*(((M - a0 - ag*G)^2)/(s1_2)) - (1/2)*log(2*pi*s1_2)
        
        A[a4] = A[a4]+A4
    }
    
    # both E and M mising
    Y=G=M=Z_age=Z_gender=NULL;A5=NULL
    a5=NULL
    a5=which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]
    
    if(length(a5)==0)    A=A
    if(length(a5)>0)
    {
        
        Y=data_matrix[a5,1]
        G=data_matrix[a5,3]
        Z_gender=data_matrix[a5,2]
        
        A5=-log(1+exp(-(((Y*pi/sqrt(3))*(b0 + bzg*Z_gender + bg*G + be*(g0 + gg*G) + (a0 + ag*G)*(be*gm + bm)))/sqrt((pi*pi/3) + (be^2)*s2_2 + s1_2* (be*gm + bm)^2))))
        
        A[a5]=A[a5]+A5
    }
    
    sum(A)
    #B=c(B,sum(A))
}







#x=seq(2,10,0.1)
#B=NULL
#for (ii in 1:length(x))
#B=c(B,f(x[ii],data_matrix))
#plot(x,B,type="l")



########################################################
########################################################
###                                                  ###
###   Calculation of log Likelihood under Null Hyp   ###
###                                                  ###
########################################################
########################################################

f_null=function(X,data_matrix)
{
    #bza=0.01;bzg=0.01
    #ag=2.4
    #gg=0.6;gm=2.3
    #s1_2=1;s2_2=1;p=0.2
    #a0=1.3;g0=1.9;b0=1
    
    bzg=X[1];
    ag=X[2];
    gg=X[3];gm=X[4];
    s1_2=X[5];
    
    s2_2=X[6];
    
    p=X[7];a0=X[8];g0=X[9];b0=X[10]
    
    if(s1_2==0) s1_2=s1_2+0.00001
    if(s2_2==0) s2_2=s2_2+0.00001
    
    if(p==0) p=p+0.00001
    if(p==1) p=p-0.00001
    
    A=NULL
    
    G=NULL; A1=NULL
    
    G=data_matrix[,3]
    
    
    A1= G*log(p) + (2-G)*log(1-p)
    A=A1
    
    # only M missing
    Y=E=G=M=Z_gender=NULL;A2=NULL
    a2=NULL
    a2=which(is.na(data_matrix[,4]=="NA")==TRUE)[which(is.na( match( which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]))==TRUE)]
    
    if(length(a2)==0)    A=A
    if(length(a2)>0)
    {
        Y=data_matrix[a2,1]
        G=data_matrix[a2,3]
        E=data_matrix[a2,5]
        Z_gender=data_matrix[a2,2]
        
        
        A2=-log(1+exp(-((Y*((s2_2 + gm^2+s1_2)*(b0 + bzg*Z_gender)))/(sqrt(s2_2 + gm^2*s1_2)*(s2_2 + gm^2*s1_2)))))- (((1/2)*(E-(g0 + gg*G + gm*(a0 + ag*G)))^2)/(s2_2+(s1_2)*(gm^2)))-(1/2)*log(2*pi*(s2_2+(gm^2)*(s1_2)))
        
        
        A[a2] = A[a2]+A2
    }
    
    
    
    #- (((1/2)*(E-(g0 + gg*G + gm*(a0 + ag*G)))^2)/(s2_2 +(gm^2*s1_2)/2))-log(2*pi*sqrt(s1_2))-(1/2)*log(2*s2_2 + gm^2 * s1_2)
    
    
    #n1=length(which(is.na(data_matrix[,7])=="FALSE"))
    
    #A2=- (log(1+exp(-(Y*(b0+bg*G + bm*M +be*E+ bza*Z_age + bzg*Z_gender)))))- (1/(2*s2_2))*((E - g0 -gg*G - gm*M)^2)- (1/2)*log(2*pi*s2_2)-(1/2)*log(2*pi*s1_2)-(1/(2*s1_2))*((M-a0-ag*G)^2)
    
    
    # only E missing
    Y=G=M=Z_gender=NULL;A3=NULL
    a3=NULL
    a3=which(is.na(data_matrix[,5]=="NA")==TRUE)[which(is.na( match( which(is.na(data_matrix[,5]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]))==TRUE)]
    
    if(length(a3)==0)    A=A
    if(length(a3)>0)
    {
        Y=data_matrix[a3,1]
        G=data_matrix[a3,3]
        M=data_matrix[a3,4]
        Z_gender=data_matrix[a3,2]
        
        A3=-log(1+exp(-((Y*(pi/sqrt(3))*(b0  + bzg*Z_gender ))/sqrt((pi*pi/3) )))) - (1/2)*(((M - a0 - ag*G)^2)/s1_2) - (1/2)*log(2*pi*s1_2)
        
        A[a3] = A[a3]+A3
        
    }
    
    
    # none of E or M missing
    Y=E=G=M=Z_gender=NULL;A4=NULL
    a4=NULL
    a4=-sort(unique(c(which(is.na(data_matrix[,5])=="TRUE"),which(is.na(data_matrix[,4])=="TRUE"))))
    if(length(which(is.na(data_matrix[,5])=="TRUE"))==0 && length(which(is.na(data_matrix[,4])=="TRUE"))==0)
    {
        for(i in 1:length(data_matrix[,5]))
        a4=c(a4,i)
    }
    
    
    
    if(length(a4)==0)    A=A
    if(length(a4)>0)
    {
        Y=data_matrix[a4,1]
        G=data_matrix[a4,3]
        M=data_matrix[a4,4]
        E=data_matrix[a4,5]
        Z_gender=data_matrix[a4,2]
        
        A4=  -log(1+exp(-(Y*(b0 + bzg*Z_gender )))) - (1/2)*(((E- g0 - gm*G - gm*M)^2)/s2_2) -(1/2)*log(2*pi*s2_2) - (1/2)*(((M - a0 - ag*G)^2)/(s1_2)) - (1/2)*log(2*pi*s1_2)
        
        A[a4] = A[a4]+A4
        
    }
    
    
    
    # both E and M mising
    Y=G=M=Z_gender=NULL;A5=NULL
    a5=NULL
    a5=which(is.na(data_matrix[,5]=="NA")==TRUE)[as.numeric(na.omit(match(which(is.na(data_matrix[,4]=="NA")==TRUE),which(is.na(data_matrix[,5]=="NA")==TRUE))))]
    
    
    if(length(a5)==0)    A=A
    if(length(a5)>0)
    {
        Y=data_matrix[a5,1]
        G=data_matrix[a5,3]
        Z_gender=data_matrix[a5,2]
        A5= -log(1+exp(-Y*(b0 + bzg*Z_gender)))
        
        A[a5] = A[a5]+A5
    }
    
    
    sum(A)
    
}

#####################################################################################
#####################################################################################
###                                                                               ###
###    Maximisation of log likelihood under Alt using optim() and EM algorithm    ###
###                                                                               ###
###                           H1: bg!=0,bm!=0,be!=0                               ###
###                                                                               ###
#####################################################################################
#####################################################################################



#mat1=cbind(mat,cov,geno,methyl,expr)
#cat(sum(pchisq(D,3,lower.tail=F)<0.05)/repli,"\n")



# Data matrix generated under Alternative
f1=function(data_matrix)
{
    #d=data_matrix[which(is.na(data_matrix[,7])=="FALSE"),] # Data matrix ignoring the individuals with missing values
    
    d=data_matrix[-sort(unique(c(which(is.na(data_matrix[,4])=="TRUE"),which(is.na(data_matrix[,5])=="TRUE")))),]
    
    
    # Initial values
    p=(2*length(which(d[,3]==2))+length(which(d[,3]==1)))/(2*dim(d)[1])
    r=cor(d[,4],d[,5])
    
    # generates alpha and sigma1^2
    meth=NULL
    meth=lm(d[,4]~d[,3])
    alphas=meth$coefficient[-1]
    s1_2=sum((meth$residuals)^2)/meth$df
    
    # generates gamma and sigma2^2
    expr=NULL
    expr=lm(d[,5]~d[,3]+d[,4])
    gammas=expr$coefficient[-1]
    s2_2=sum((expr$residuals)^2)/expr$df
    
    # generates beta
    y=glm(((d[,1]+1)/2)~d[,2]+d[,3]+d[,4]+d[,5],family=binomial)
    betas=y$coefficient[-1]
    if (max(as.numeric(betas))<10^10)
    {
        #par_actual=c(beta.age,beta.gender,beta.g,beta.m,beta.e,alpha.g,gamma.g,gamma.m,sig.m^2,sig.e^2,MAF,alpha0,gamma0,beta0)
        par_initial=c(as.numeric(betas),as.numeric(alphas),as.numeric(gammas),s1_2,s2_2,p,as.numeric(meth$coefficient[1]),as.numeric(expr$coefficient[1]),as.numeric(y$coefficient[1]))
        
        
        mat=NULL
        #mat=par_actual
        #rownames(mat)=c("bza","bzg","bg","bm","be","ag","gg","gm","s1_2","s2_2","p","a0","g0","b0")
        mat=par_initial
        
        op=optim(par_initial,f,data_matrix=data_matrix,method="L-BFGS-B",lower=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf,Inf),control = list(fnscale=-1,maxit=500))
        
        mat=cbind(mat,as.numeric(op$par))
        k=0
        while(k<100)
        {
            if(max(abs(mat[,dim(mat)[2]]-mat[,(dim(mat)[2]-1)]))>0.00001)
            {
                op=optim(mat[,dim(mat)[2]],f,data_matrix=data_matrix,method="L-BFGS-B",lower=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf,Inf),control = list(fnscale=-1,maxit=500))
                
                mat=cbind(mat,as.numeric(op$par))
                k=k+1
                
            }
            if(max(abs(mat[,dim(mat)[2]]-mat[,(dim(mat)[2]-1)]))<0.00001)
            break;
        }
        
        #return(mat[c(1:3,length(mat)[2])])
        return(mat)
        
    }
    
}



#####################################################################################
#####################################################################################
###                                                                               ###
###    Maximisation of log likelihood under Null using optim() and EM algorithm   ###
###                                                                               ###
###                                  H0: bg=bm=be=0                               ###
###                                                                               ###
#####################################################################################
#####################################################################################

# Data matrix generated under Null Hypothesis: bg=bm=be=0

f2=function(data_matrix)
{
    
    #d=data_matrix[which(is.na(data_matrix[,5])=="FALSE"),] # Data matrix ignoring the individuals with missing values
    d=data_matrix[-sort(unique(c(which(is.na(data_matrix[,4])=="TRUE"),which(is.na(data_matrix[,5])=="TRUE")))),]
    
    # Initial values
    p=(2*length(which(d[,3]==2))+length(which(d[,3]==1)))/(2*dim(d)[1])
    r=cor(d[,4],d[,5])
    
    # generates alpha and sigma1^2
    meth=lm(d[,4]~d[,3])
    alphas=meth$coefficient[-1]
    s1_2=sum((meth$residuals)^2)/meth$df
    
    # generates gamma and sigma2^2
    expr=lm(d[,5]~d[,3]+d[,4])
    gammas=expr$coefficient[-1]
    s2_2=sum((expr$residuals)^2)/expr$df
    
    # generates beta
    y=glm(((d[,1]+1)/2)~d[,2],family=binomial)
    betas=y$coefficient[-1]
    if (max(as.numeric(betas))<10^10)
    {
        #par_actual=c(beta.age,beta.gender,alpha.g,gamma.g,gamma.m,sig.m^2,sig.e^2,MAF,alpha0,gamma0,beta0)
        par_initial=as.numeric(c(as.numeric(betas)[1],as.numeric(alphas),as.numeric(gammas),s1_2,s2_2,p,as.numeric(meth$coefficient[1]),as.numeric(expr$coefficient[1]),as.numeric(y$coefficient[1])))
        
        
        
        mat=NULL
        mat=par_initial
        
        
        op=optim(par_initial,f_null,data_matrix=data_matrix,method="L-BFGS-B",lower=c(-Inf,-Inf,-Inf,-Inf,0,0,0,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf,Inf),control = list(fnscale=-1,maxit=500))
        
        
        mat=cbind(mat,as.numeric(op$par))
        
        k=0
        while(k<100)
        {
            if(max(abs(mat[,dim(mat)[2]]-mat[,(dim(mat)[2]-1)]))>0.00001)
            {
                op=optim(mat[,dim(mat)[2]],f_null,data_matrix=data_matrix,method="L-BFGS-B",lower=c(-Inf,-Inf,-Inf,-Inf,0,0,0,-Inf,-Inf,-Inf), upper=c(Inf,Inf,Inf,Inf,Inf,Inf,1,Inf,Inf,Inf),control = list(fnscale=-1,maxit=500))
                
                #op=optim(mat[,dim(mat)[2]],f,data_matrix=data_matrix,method="SANN",control = list(fnscale=-1,maxit=500))
                
                mat=cbind(mat,as.numeric(op$par))
                #cat(dim(mat)[2],"\t")
                #cat(max(abs(mat[,dim(mat)[2]]-mat[,(dim(mat)[2]-1)])),"\n")
                k=k+1
                
            }
            if(max(abs(mat[,dim(mat)[2]]-mat[,(dim(mat)[2]-1)]))<0.00001)
            break;
        }
        
        #return(mat[,c(1:3,dim(mat)[2])])
        return(mat)
        
    }
    
    
    
}


############################  Likelihood analysis code ends    ####################
exprn=as.numeric(bb[1,-(1:5)])
geno.j=as.numeric(bb[2,-(1:5)])
meth.i=as.numeric(bb[3,-(1:5)])

matid=read.table("mat.id.txt",sep="\t")
mat.id=as.character(matid[,1])
covariates_file=read.table("covariates_file.txt",sep="\t",header=T)

cov=covariates_file[c(which(covariates_file[,6]==0),which(covariates_file[,6]==1)),c(6,5,4)]
cov[which(cov[,1]==0),1]=-1


###masterfunc=function()
###{
    D=data.frame(matrix(nrow=1,ncol=8))
    mat1=NULL
    mat1=cbind(mat.id,cov,geno.j,meth.i,exprn)
    mat1[,1]=as.numeric(unlist(strsplit(as.character(mat1[,1]),split="Subject_ID_"))[seq(2,(2*(dim(mat1)[1])),by=2)])
    mat1[,4]=as.numeric(mat1[,4]) # 1 denotes female; 2 denotes male
    mat1[,5]=mat1[,5]-1
    colnames(mat1)[1:2]=paste(c("Subject_ID","Status"),sep="")
    colnames(mat1)[5:6]=paste(c("geno","methyl"),sep="")
    data_matrix1=NULL
    data_matrix1=mat1[,-c(1,3)]
    check_geno=data_matrix1[-sort(unique(c(which(is.na(data_matrix1[,4])=="TRUE"),which(is.na(data_matrix1[,5])=="TRUE")))),3]
    if(length(unique(check_geno))==1) cat()
    if(length(unique(check_geno))>1)
    {
        x=x_null=NULL
        func=func_null=NULL
        func=tryCatch({suppressWarnings(f1(data_matrix1))},error=function(e){})
        if (length(func)==0) cat()
        if (length(func)>0)
        {
            x=as.vector(func[,dim(func)[2]])
            alt_lik=f(x,data_matrix1)

            func_null=tryCatch({suppressWarnings(f2(data_matrix1))},error=function(e){})
            if (length(func_null)==0) cat()
            if (length(func_null)>0)
            {
                x_null=as.vector(func_null[,dim(func_null)[2]])
                null_lik=f_null(x_null,data_matrix1)
                
                #mi=which(meth.case[,23]==Gene.name)[meth.i]
            
                D=cbind(2*(alt_lik-null_lik),pchisq(2*(alt_lik-null_lik),3,lower.tail=F),genesymbol,paste(as.character(bb[3,5])),bb[2,1:4])
                colnames(D)=paste(c("stat","pval","GENESYMBOL","CpG","Name","SNP","Chr","MapInfo"))
                
           
                #D=cbind(2*(alt_lik-null_lik),pchisq(2*(alt_lik-null_lik),3,lower.tail=F),paste(meth.case[mi,23]),paste(meth.case[mi,1]),meth.case[mi,25],paste(case.geno.file.gene[geno.j,1]),paste(case.geno.file.gene[geno.j,2]),case.geno.file.gene[geno.j,3],case.geno.file.gene[geno.j,4])
                ###return(D)
                
            }
            
        }
    }
    
###}

###bb=masterfunc(i,j)

### End of main code
#colnames(D)=paste(c("stat","pval","GENESYMBOL","CpG","Name","SNP","Chr","MapInfo"))
write.table(D,fout,row.names=F,col.names=T,quote=F,sep="\t")

#  j: genotype;   i: methylation


