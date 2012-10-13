% File man/NewInput.MMC.Rd
\name{NewInput.MMC}
\alias{NewInput.MMC}
\title{Define the inputs of a new M/M/c queueing model}
\description{
  Define the inputs of a new M/M/c queueing model
}

\usage{
  NewInput.MMC(lambda=0, mu=0, c=1, n=0, method=0)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{c}{number of servers}
  \item{n}{number of clients in the system. Put n=0 for a idle probability. With n=-1, no probabilities are computed}
  \item{method}{method of computation of the probabilities of n number of clients in the system. With method=0, the exact results are calculated using the formal definition. With method=1, aproximate results are calculated using Stirling aproximation of factorials and logaritms.}
}

\details{Define the inputs of a new M/M/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMC}}
}


\examples{
## See example 10.9 in reference [Sixto2004] for more details.
## create input parameters
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0, method=0)
}
\keyword{M/M/c}

