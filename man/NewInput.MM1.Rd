% File man/NewInput.MM1.Rd
\name{NewInput.MM1}
\alias{NewInput.MM1}
\title{Define the inputs of a new M/M/1 queueing model}
\description{
  Define the inputs of a new M/M/1 queueing model
}

\usage{
NewInput.MM1(lambda=0, mu=0, n=0)
}

\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{n}{number of customers in the system. Put n=0 for a idle probability. With n=-1, no probabilities are computed}
}

\details{Define the inputs of a new M/M/1 queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{CheckInput.i_MM1}}
}


\examples{
## See example 10.3 in reference [Sixto2004] for more details.
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)
}

\keyword{M/M/1}
