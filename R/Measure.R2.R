Measure.R2 <-
function(biloci){
M.r2		 =	NA

if (any(is.na(biloci)==TRUE)){

	ligne		=	na.action(na.omit(biloci))
	biloci		=	biloci[-ligne,]

}
SIG=var(biloci)

ifelse((SIG[1,1]<0.0000001) | (SIG[2,2]<0.0000001),M.r2<-0, M.r2<-(SIG[1,2])^2/(SIG[1,1]*SIG[2,2]) )
as.numeric(M.r2)

}

