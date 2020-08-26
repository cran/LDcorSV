\name{Measure.R2S}
\alias{Measure.R2S}
\title{r^2_S measure}

\description{This function estimates the novel measure of linkage disequilibrium which is corrected by the structure of the sample.}

\usage{Measure.R2S(biloci, struc, na.presence=TRUE)}

\arguments{

  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes)

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Row names correspond to the  ID of individuals.

Column names correspond to the ID of markers.}

  \item{struc}{Numeric matrix (N x (P-1)), where N is the number of genotypes (or haplotypes) and P the number of sub-populations.

Matrix values are the probabilities for each genotypes (or haplotypes) to belong to each sub-populations.

Row names must correspond to the ID of individuals and must be ranged as in the biloci matrix.

Column names correspond to the ID of sub-populations.

The matrix must be inversible, if the structure is with P sub-populations, only P-1 columns are expected.

No missing value.}

  \item{na.presence}{Boolean indicating the presence of missing values in data.
  
  If na.presence=FALSE (no missing data), computation of \eqn{r^2_V} and \eqn{r^2_{VS}} is largely optimized.
  
  By default, na.presence=TRUE.}

}

\value{The returned value is the estimated value of the measure of linkage disequilibrium  corrected by  the structure of the sample or NA if less than 5 individuals have non-missing data at both loci.}

\references{Mangin, B., Siberchicot, A., Nicolas, S., Doligez, A., This, P., Cierco-Ayrolles, C. (2012). Novel measures of linkage disequilibrium that correct the bias due to population structure and relatedness. Heredity, 108 (3), 285-291. DOI: 10.1038/hdy.2011.73}

\author{David Desrousseaux, Florian Sandron, Aurélie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin}

\examples{
data(data.test)
Geno <- data.test[[1]]
S.2POP <- data.test[[3]]
Measure.R2S(Geno, S.2POP)
}