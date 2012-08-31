\name{Measure.R2V}

\alias{Measure.R2V}

\title{r^2_V measure}

\description{This function estimates the novel measure of linkage disequilibrium which is corrected by the relatedness of genotyped individuals.}

\usage{Measure.R2V(biloci, V)}

\arguments{

  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes).

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Row names correspond to the  ID of individuals.

Column names correspond to the  ID of markers.}

  \item{V}{Numeric matrix (N x N), where N is the number of genotypes (or haplotypes).
	
	Matrix values are coefficients of genetic covariance for each pair of individuals.

	Rows and columns names must correspond to the ID of individuals and must be ranged in the same order as in the biloci matrix.

No missing value.}

	}

\value{The returned value is the estimated value of the measure of linkage disequilibrium corrected by the relatedness of genotyped individuals or NA if less than 5 individuals have non-missing data at both loci.}


\references{Mangin et al, 2012,  Novel measures of linkage disequilibrium that correct for the bias due to population structure and relatedness.}
\author{David Desrousseaux, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin }





