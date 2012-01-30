\name{LD.Measures}
\alias{LD.Measures}

\title{LD Measures}

\description{This function estimates for a set of loci:

- the usual measure of linkage disequilibrium (\eqn{r^2}).

- the measure of linkage disequilibrium corrected by the structure of the sample (\eqn{r^2_S}).

- the measure of linkage disequilibrium corrected by the relatedness of genotyped individuals (\eqn{r^2_V}).

- the measure of linkage disequilibrium corrected  by both, the relatedness of genotyped individuals and the structure of the sample  (\eqn{r^2_{VS}}).

This function gives extra informations on the studied loci.}

\usage{LD.Measures(donnees, V = NA, S = NA, data = "G", supinfo = FALSE)}

\arguments{
  \item{donnees}{Numeric matrix (N x M), where N is the number of genotypes (or haplotypes) and M is the number of markers.

Matrix values are the allelic doses:

	- (0,1,2) for genotypes.

	- (0,1) for haplotypes.


Row names correspond to the ID of individuals.

Column names correspond to the ID  of markers.

Missing values are allowed.}

  \item{V}{Numeric matrix (N x N), where N is the number of genotypes (or haplotypes).
	
	Matrix values are coefficients of genetic variance-covariance between every pair of individuals.

	Row and column names must correspond to the ID of individuals.

No missing value.}

  \item{S}{Numeric matrix (N x (1-P)), where N is the number of genotypes (or haplotypes) and P the number of sub-populations.

Matrix values are the probabilities (between 0 and 1) for each genotype (or haplotype) to belong to each sub-populations.

Row names must correspond to the ID of individuals.

Column names correspond to the ID of sub-populations.

The matrix must be inversible, if the structure is with P sub-populations, only P-1 columns are expected.

No missing value.}

  \item{data}{Value equal to "G" or "H" depending on the type of data (Genotype or Haplotype).

	Default value is "G".}


  \item{supinfo}{Boolean indicating whether you wish to get information about the loci.

	If Supinfo=TRUE, for each locus, the Minor Allelic Frequency (MAF), the frequency of heterozygous genotypes (only if the data are genotypes) and the missing value frequency are computed.

 By default, supinfo=FALSE.}

	}

\value{The returned value is a dataframe of size (M(M-1))/2 rows and C columns, where M is the number of markers  and  C is a number between 3 and 12 depending on options chosen by user.

  The first three columns contain respectively the name of the first marker, the name of the second marker and the estimated value of the usual measure of linkage disequilibrium (\eqn{r^2}) between these two markers.

If  only V is different from  NA, the fourth column contains the estimated value of the measure of linkage disequilibrium corrected by the relatedness of genotyped individuals (\eqn{r^2_V}).

If  only S is different from  NA, the fourth column contains the estimated value of the measure of linkage disequilibrium corrected by relatedness corrected by the structure of the sample (\eqn{r^2_S}).

If V and S are simultaneously different from NA, the fourth, fifth and sixth columns respectively contain the estimated values of \eqn{r^2_V}, \eqn{r^2_S} and \eqn{r^2_{VS}} (\eqn{r^2} corrected  by both the relatedness of genotyped individuals and the structure of the sample).

If Supinfo=TRUE, then the last six columns contain information for both loci : the MAF, the frequency of heterozygous genotype (NA if haplotype data) and the missing value frequency.}

\references{Mangin et al, 2011,  Novel Measures of linkage Disequilibrium that correct for the bias due to population structure and relatedness, Heredity.}

\author{David Desrousseaux, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin }




\examples{
data(data.test)
Geno<-data.test[[1]]
Geno.test<-Geno[,1:5]
V.WAIS<-data.test[[2]]
S.2POP<-data.test[[3]]
LD<-LD.Measures(Geno.test,V=V.WAIS,S=S.2POP,data ="G",supinfo=TRUE)
LD
}


