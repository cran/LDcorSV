\name{Inv.proj.matrix.sdp}
\alias{Inv.proj.matrix.sdp}

\title{Inv.proj.matrix.sdp}

\description{This function computes the Moore-Penrose pseudo-inverse of a symetric matrix. A single value decomposition is performed, the non positive eigen values are set to zero, then the pseudo-inverse is computed.}
\usage{Inv.proj.matrix.sdp(matrix)}


\arguments{
  \item{matrix}{symmetric matrix}
	  }

\value{The returned value is the pseudo-inverse matrix.}


\author{ David Desrousseaux, Aurelie Siberchicot, Christine Cierco-Ayrolles and Brigitte Mangin }





\examples{
data(data.test)
V.WAIS<-data.test[[2]]
Inv.V.WAIS<-Inv.proj.matrix.sdp(V.WAIS)
Inv.V.WAIS
}



