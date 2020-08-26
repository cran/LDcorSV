\name{Measure.R2}
\alias{Measure.R2}
\title{r^2 measure}

\description{This function estimates the usual measure of linkage disequilibrium (\eqn{r^2})  between two loci.}

\usage{Measure.R2(biloci, na.presence = TRUE)}

\arguments{

  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes)

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Rows names correspond to the ID of individuals.

Columns names correspond to the ID of markers.}

  \item{na.presence}{Boolean indicating the presence of missing values in data.
  
  If na.presence=FALSE (no missing data), computation of \eqn{r^2_V} and \eqn{r^2_{VS}} is largely optimized.
  
  By default, na.presence=TRUE.}
}

\value{
    The returned value is the estimated value of the usual measure of linkage disequilibrium (\eqn{r^2}) 
    or NA if less than 5 individuals have non-missing data at both loci
}

\references{Hill, W.G, Robertson, A. (1968). Linkage diseqilibrium in finite populations. Theoretical and Applied Genetics, 38, 226-231. DOI: 10.1007/BF01245622}

\author{David Desrousseaux, Florian Sandron, Aurélie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin}


\examples{
data(data.test)
Geno <- data.test[[1]]
Measure.R2(Geno)
}