% File man/NewInput.MM1K.Rd
\name{NewInput.MM1K}
\alias{NewInput.MM1K}
\title{Define the inputs of a new M/M/1/K queueing model}
\description{
  Define the inputs of a new M/M/1/K queueing model
}

\usage{
  \method{NewInput}{MM1K}(lambda=0, mu=0, k=1)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{k}{system capacity}
}

\details{Define the inputs of a new M/M/1/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1K}}
}

\examples{
## See example 10.7 in reference [Sixto2004] for more details.
## create input parameters
i_mm1k <- NewInput.MM1K(lambda=5, mu=5.714, k=15)

}
\keyword{M/M/1/K}

