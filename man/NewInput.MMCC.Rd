% File man/NewInput.MMCC.Rd
\name{NewInput.MMCC}
\alias{NewInput.MMCC}
\title{Define the inputs of a new M/M/c/c queueing model}
\description{
  Define the inputs of a new M/M/c/c queueing model
}

\usage{
  \method{NewInput}{MMCC}(lambda=0, mu=0, c=1)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{c}{number of servers}
}

\details{Define the inputs of a new M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMCC}}
}



\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

}
\keyword{M/M/c/c}

