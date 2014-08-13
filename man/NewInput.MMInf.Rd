% File man/NewInput.MMInf.Rd
\name{NewInput.MMInf}
\alias{NewInput.MMInf}
\title{Define the inputs of a new M/M/Infinite queueing model}
\description{
  Define the inputs of a new M/M/Infinite queueing model
}
\usage{
NewInput.MMInf(lambda=0, mu=0, n=0)
}

\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{n}{number of customers in the system. Put n=0 for a standard model}
}

\details{Define the inputs of a new M/M/Infinite queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MMInf}}
}


\examples{
## create input parameters
i_mminf <- NewInput.MMInf(lambda=0.25, mu=4, n=0)

}
\keyword{M/M/Infinite}
