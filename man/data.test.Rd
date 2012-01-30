\name{data.test}
\alias{data.test}
\docType{data}
\title{data.test}

\description{"data.test" is a list of 3 elements: 

- \code{Geno}: allelic doses of 20 markers on a chromosome of 91 Vitis vinifera plants.

- \code{V.WAIS}: kinship matrix of 91 plants of Vitis vinifera.

- \code{S.2POP}: structure population matrix of 91 plants of Vitis vinifera in two sub-populations.}

\usage{data(data.test)}


\format{

A list containing the following components:

- \code{Geno}: matrix (91 x 20) of numerical values

- \code{V.WAIS}: matrix (91 x 91) of numerical values

- \code{S.2POP}: matrix (91 x 1) of  numerical values


}
\details{}


\examples{

data(data.test)

#Allelic doses of 20 markers on a chromosome of 91  Vitis vinifera plants
Geno=data.test[[1]]
Geno

#Kinship matrix of 91 plants of Vitis vinifera
V.WAIS=data.test[[2]]
V.WAIS

#Structure population matrix of 91 plants of Vitis vinifera 
# in two sub-populations
S.2POP=data.test[[3]]
S.2POP

}
\keyword{datasets}
