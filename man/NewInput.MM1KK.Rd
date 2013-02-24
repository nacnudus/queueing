% File man/NewInput.MM1KK.Rd
\name{NewInput.MM1KK}
\alias{NewInput.MM1KK}
\title{Define the inputs of a new M/M/1/K/K queueing model}
\description{
  Define the inputs of a new M/M/1/K/K queueing model
}

\usage{
  NewInput.MM1KK(lambda=0, mu=0, k=1, method=0)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{k}{system capacity}
  \item{method}{method of computation of the probabilities of k (system capacity) clients down. With method=0, the exact results are calculated using the formal definition. With method=1, aproximate results are calculated using Stirling aproximation of factorials and logaritms. With method=2, Jain's Method [Jain2007], pag. 26 is used}
}

\details{Define the inputs of a new M/M/1/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
\cr
  [Jain2007] Joti Lal Jain, Sri Gopal Mohanty, Walter Bohm (2007).\cr
\emph{A course on Queueing Models}.\cr
Chapman-Hall.
}

\seealso{
  \code{\link{CheckInput.i_MM1KK}}
}

\examples{
## See example 10.13 in reference [Sixto2004] for more details.
## create input parameters
i_mm1kk <- NewInput.MM1KK(lambda=0.25, mu=4, k=2, method=0)

}
\keyword{M/M/1/K/K}
