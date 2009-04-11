% File man/NewInput.MMCK.Rd
\name{NewInput.MMCK}
\alias{NewInput.MMCK}
\title{Define the inputs of a new M/M/c/K queueing model}
\description{
  Define the inputs of a new M/M/c/K queueing model
}

\usage{
  \method{NewInput}{MMCK}(lambda=0, mu=0, c=1, k=1)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{c}{number of servers}
  \item{k}{system capacity}
}

\details{Define the inputs of a new M/M/c/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCK}}
}


\examples{
## See example 10.11 in reference [Sixto2004] for more details.
## create input parameters
i_mmck <- NewInput.MMCK(lambda=8, mu=4, c=5, k=12)

}
\keyword{M/M/c/K}

