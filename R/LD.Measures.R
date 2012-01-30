LD.Measures <-
function(donnees,V=NA,S=NA,data="G",supinfo=FALSE){


library(MASS)


loc1<-vector("character",0)
loc2<-vector("character",0)
r2<-numeric(0)
M.r2<-NA	
r2v<-numeric(0)
M.r2v<-NA
r2s<-numeric(0)
M.r2s<-NA
r2vs<-numeric(0)
M.r2vs<-NA
MAF.loc1<-numeric(0)
heterofreq.loc1<-numeric(0)
NAfreq.loc1<-numeric(0)
MAF.loc2<-numeric(0)
heterofreq.loc2<-numeric(0)
NAfreq.loc2<-numeric(0)

is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

if(is.null(rownames(donnees))==TRUE | is.null(colnames(donnees))==TRUE) 
{return(print("ERROR: genotype data must have row and column names"))}
mark<-colnames(donnees)


if(all(is.wholenumber(donnees[!is.na(donnees)])==FALSE))
{return(print("ERROR: haplotype or genotype data must be a numeric matrix with integer values"))}

if( data=="G"){if(max(donnees[!is.na(donnees)])>2 | min( donnees[!is.na(donnees)] )<0)
{return(print("ERROR: genotype data must be a numeric matrix with integer values between 0 and 2"))}}


if( data=="H"){if(max(donnees[!is.na(donnees)])>1 | min(donnees[!is.na(donnees)])<0)
{return(print("ERROR: haplotype must be a numeric matrix with integer values between 0 and 1"))}}

if(length(unique(is.na(S)))==2)
{return(print("ERROR: missing values are not allowed in the structure matrix"))}

if(length(unique(is.na(V)))==2)
{return(print("ERROR: missing values are not allowed in the genetic variance-covariance matrix"))}

S.yes=all(is.na(S))
V.yes=all(is.na(V))

if (S.yes==FALSE) {
if( max(S)>1 | min(S)<0)
	{return(print("ERROR: the structure matrix must be a numeric matrix with values between 0 and 1"))}
if( is.null(colnames(S))==TRUE |  is.null(rownames(S))==TRUE) 
{return(print("ERROR:  structure matrix  must have row and column names"))}
                   }
if(V.yes==FALSE) {
if( is.null(colnames(V))==TRUE |  is.null(rownames(V))==TRUE) 
{return(print("ERROR:  genotypic variance-covariance matrix  must have row and column names"))}
                   }

if (S.yes==FALSE & V.yes==FALSE )
{
ID=rownames(donnees)
ID.S=rownames(S)
ID.Vr=rownames(V)
ID.Vc=colnames(V)
     if(length(setdiff(ID,ID.S))!=0)
	{return(print("ERROR: some ID genotype are not found in ID structure "))}
     if(length(setdiff(ID,ID.Vr))!=0 | length(setdiff(ID,ID.Vc))!=0 )
	{return(print("ERROR: some ID genotype are not found in ID genotypic variance-covariance"))}



			alignS	          =	merge(donnees,S,by.x=0,by.y=0,all.x=TRUE,sort=FALSE)
			donnees            =     alignS[,2:(ncol(donnees)+1)]
			rownames(donnees)  =     alignS[,1]
			
			S           	  =     as.data.frame(alignS[,(ncol(donnees)+2):ncol(alignS)])
			rownames(S)       =     alignS[,1]
			colnames(S)       =	colnames(alignS)[(ncol(donnees)+2):ncol(alignS)]


			alignV		=	merge(V,donnees,by.x=0,by.y=0,sort=FALSE)
			donnees		=	alignV[,(ncol(V)+2):ncol(alignV)]
			rownames(donnees)=	alignV[,1]

			rownames(alignV)=	alignV[,1]
			V		=	alignV[,-1]
			V		=	t(V)
			V		=	merge(V,donnees,by.x=0,by.y=0,sort=FALSE)
			rownames(V)	=	V[,1]
			V		=	V[,2:(nrow(V)+1)]


}
else
{
if (S.yes==FALSE)
{
ID=rownames(donnees)
ID.S=rownames(S)
     if(length(setdiff(ID,ID.S))!=0)
	{return(print("ERROR: some ID genotype are not found in ID structure "))}
 


			alignS	          =	merge(donnees,S,by.x=0,by.y=0,all.x=TRUE,sort=FALSE)
			donnees            =     alignS[,2:(ncol(donnees)+1)]
			rownames(donnees)  =     alignS[,1]
			
			S           	  =     as.data.frame(alignS[,(ncol(donnees)+2):ncol(alignS)])
			rownames(S)       =     alignS[,1]
			colnames(S)       =	colnames(alignS)[(ncol(donnees)+2):ncol(alignS)]



}

if (V.yes==FALSE)

{
ID=rownames(donnees)
ID.Vr=rownames(V)
ID.Vc=colnames(V)
     if(length(setdiff(ID,ID.Vr))!=0 | length(setdiff(ID,ID.Vc))!=0 )
	{return(print("ERROR: some ID genotype are not found in ID genotypic variance-covariance"))}
alignV		=	merge(V,donnees,by.x=0,by.y=0,sort=FALSE)
			donnees		=	alignV[,(ncol(V)+2):ncol(alignV)]
			rownames(donnees)=	alignV[,1]

			rownames(alignV)=	alignV[,1]
			V		=	alignV[,-1]
			V		=	t(V)
			V		=	merge(V,donnees,by.x=0,by.y=0,sort=FALSE)
			rownames(V)	=	V[,1]
			V		=	V[,2:(nrow(V)+1)]



}
}


if (supinfo==TRUE) {info<-apply(donnees,data=data,2,Info.Locus)}

for (i in 1:(ncol(donnees)-1)){	
	for (j in (i+1):ncol(donnees)){

li<-donnees[,i]
lj<-donnees[,j]
lij<-cbind(li,lj)
loc1<-c(loc1,mark[i])
loc2<-c(loc2,mark[j])



options(warn=-1)
M.r2<-Measure.R2(biloci=lij)
options(warn=0)			

#print(c("i=",i,"j=",j))

if (S.yes==FALSE){M.r2s<-Measure.R2S(biloci=lij,struc=S)}
if (V.yes==FALSE){M.r2v<-Measure.R2V(biloci=lij,V=V)}
if (S.yes==FALSE & V.yes==FALSE){M.r2vs<-Measure.R2VS(biloci=lij,V=V,struc=S)}
			

r2<-c(r2,M.r2)
r2s<-c(r2s,M.r2s)
r2v<-c(r2v,M.r2v)
r2vs<-c(r2vs,M.r2vs)		



if (supinfo==TRUE){M.MAF.loc1<-NA
		   M.heterofreq.loc1<-NA
		   M.NAfreq.loc1<-NA
	           M.MAF.loc2<-NA
		   M.heterofreq.loc2<-NA
		   M.NAfreq.loc2<-NA

			supinfo.li<-info[,which(colnames(info)==mark[i])]
			supinfo.lj<-info[,which(colnames(info)==mark[j])]
			M.MAF.loc1<-supinfo.li[1]
			M.heterofreq.loc1<-supinfo.li[2]
			M.NAfreq.loc1<-supinfo.li[3]
			M.MAF.loc2<-supinfo.lj[1]
			M.heterofreq.loc2<-supinfo.lj[2]
			M.NAfreq.loc2<-supinfo.lj[3]

		

				MAF.loc1<-c(MAF.loc1,M.MAF.loc1)
				heterofreq.loc1<-c(heterofreq.loc1,M.heterofreq.loc1)
				NAfreq.loc1<-c(NAfreq.loc1,M.NAfreq.loc1)
				MAF.loc2<-c(MAF.loc2,M.MAF.loc2)
				heterofreq.loc2<-c(heterofreq.loc2,M.heterofreq.loc2)
				NAfreq.loc2<-c(NAfreq.loc2,M.NAfreq.loc2)}
				

			}

					}


if (S.yes==FALSE & V.yes==FALSE){result<-data.frame(loc1,loc2,r2,r2v,r2s,r2vs)}

else{if (S.yes==FALSE){result<-data.frame(loc1,loc2,r2,r2s)}
	else{if (V.yes==FALSE){result<-data.frame(loc1,loc2,r2,r2v)}
		else{result<-data.frame(loc1,loc2,r2)}}}

if (supinfo==TRUE){Info<-data.frame(MAF.loc1,heterofreq.loc1,NAfreq.loc1,MAF.loc2, heterofreq.loc2,NAfreq.loc2)
		   result<-c(result,Info)}

result<-as.data.frame(result)
result
}
