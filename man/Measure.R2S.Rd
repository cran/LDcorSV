\name{Measure.R2S}
\alias{Measure.R2S}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{r^2_S measure}

\description{This function estimates the novel measure of linkage disequilibrium which is corrected by the structure of the sample.}

\usage{Measure.R2S(biloci, struc)}

\arguments{
  \item{biloci}{Numeric matrix (N x 2), where N is the number of genotypes (or haplotypes)

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.

Row names correspond to the  ID of individuals.

Column names correspond to the  ID of markers.}


  \item{struc}{Numeric matrix (N x (P-1)), where N is the number of genotypes (or haplotypes) and P the number of sub-populations.

Matrix values are the probabilities for each genotypes (or haplotypes) to belong to each sub-populations.

Row names must correspond to the ID of individuals and must be ranged as in the biloci matrix.

Column names correspond to the ID of sub-populations.

The matrix must be inversible, if the structure is with P sub-populations, only P-1 columns are expected.

No missing value.}
	}


\value{The returned value is the estimated value of the measure of linkage disequilibrium  corrected by  the structure of the sample.}

\references{Mangin et al, 2011,  Novel measures of linkage disequilibrium that correct for the bias due to population structure and relatedness, Heredity.}
\author{David Desousseaux, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin }



\seealso{}




