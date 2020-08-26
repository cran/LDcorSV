\name{Measure.R2V}
\alias{Measure.R2V}
\title{r^2_V measure}

\description{This function estimates the novel measure of linkage disequilibrium which is corrected by the relatedness of genotyped individuals.}

\usage{Measure.R2V(biloci, V, na.presence=TRUE, V_inv=NULL)}

\arguments{

  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes).

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Row names correspond to the ID of individuals.

Column names correspond to the ID of markers.}

  \item{V}{Numeric matrix (N x N), where N is the number of genotypes (or haplotypes).
	
	Matrix values are coefficients of genetic covariance for each pair of individuals.

	Rows and columns names must correspond to the ID of individuals and must be ranged in the same order as in the biloci matrix.

No missing value.}

 \item{na.presence}{Boolean indicating the presence of missing values in data.
  
  If na.presence=FALSE (no missing data), computation of \eqn{r^2_V} and \eqn{r^2_{VS}} is largely optimized.
  
  By default, na.presence=TRUE.}
  
  \item{V_inv}{Should stay NULL.}
}

\value{The returned value is the estimated value of the measure of linkage disequilibrium corrected by the relatedness of genotyped individuals or NA if less than 5 individuals have non-missing data at both loci.}

\references{Mangin, B., Siberchicot, A., Nicolas, S., Doligez, A., This, P., Cierco-Ayrolles, C. (2012). Novel measures of linkage disequilibrium that correct the bias due to population structure and relatedness. Heredity, 108 (3), 285-291. DOI: 10.1038/hdy.2011.73}

\author{David Desrousseaux, Florian Sandron, Aurélie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin}

\examples{
data(data.test)
Geno <- data.test[[1]]
V.WAIS <- data.test[[2]]
Measure.R2V(Geno, V.WAIS)
}