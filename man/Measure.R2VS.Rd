\name{Measure.R2VS}
\alias{Measure.R2VS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{r^2_VS measure}
\description{This function estimates the novel measure of linkage disequilibrium which is corrected by  both the relatedness of genotyped individuals and  the structure of the sample.

}
\usage{
Measure.R2VS(biloci, V, struc, na.presence=TRUE, V_inv=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes)

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Row names correspond to the  ID of individuals.

Column names correspond to the  ID of markers.

}
  \item{V}{Numeric matrix (N x N), where N is the number of genotypes (or haplotypes).

	Matrix values are coefficients of genetic variance-covariance for every pair of individuals.
	Row and column names must correspond to the ID of individuals and must be ranged as in the biloci matrix.

No missing value.
}
  \item{struc}{Numeric matrix (N x (P-1)), where N is the number of genotypes (or haplotypes) and P the number of sub-populations.

Matrix values are the probabilities (between 0 and 1)
 for each genotypes (or haplotypes) to belong to each sub-populations.

Row names must correspond to the ID of individuals and must be ranged as in the biloci matrix.

Column names correspond to the ID of sub-populations.

The matrix must be inversible, if the structure is with P sub-populations, only P-1 columns are expected.

No missing value.}


 \item{na.presence}{Boolean indicating the presence of missing values in data.
  
  If na.presence=FALSE (no missing data), computation of \eqn{r^2_V} and \eqn{r^2_{VS}} is largely optimized.
  
  By default, na.presence=TRUE.}
  
  \item{V_inv}{Should stay NULL}
}

\value{The returned value is the estimated value of the linkage disequilibrium measure corrected by both the relatedness of genotyped individuals and  the structure of the sample or NA if less than 5 individuals have non-missing data at both loci.

}
\references{Mangin et al, 2012,  Novel measures of linkage disequilibrium that correct for the bias due to population structure and relatedness, Heredity.}

\author{David Desrousseaux, Florian Sandron, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin }






% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

