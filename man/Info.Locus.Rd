\name{Info.Locus}
\alias{Info.Locus}



\title{Information on loci}

\description{For a locus, this function computes the minor allelic frequency, the frequency of heterozygous genotypes and the missing value frequency.}

\usage{Info.Locus(locus,data="G")}

\arguments{

	\item{locus}{Numeric vector of allelic doses.}

	\item{data}{Value equal to "G" or "H" depending on the type of data (Genotype or Haplotype).
	Default value is "G".}
	 
	 }

\value{The returned value is a numeric vector of three values which are respectively the minor allelic frequency, the  frequency of heterozygous genotypes (NA if haplotype data) and the missing value frequency.}


\author{David Desrousseaux, Florian Sandron, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin}


\examples{
data(data.test)
Geno<-data.test[[1]]
info<-apply(Geno,2,Info.Locus)
info
}



