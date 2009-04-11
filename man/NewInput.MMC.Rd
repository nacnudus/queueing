% File man/NewInput.MMC.Rd
\name{NewInput.MMC}
\alias{NewInput.MMC}
\title{Define the inputs of a new M/M/c queueing model}
\description{
  Define the inputs of a new M/M/c queueing model
}

\usage{
  \method{NewInput}{MMC}(lambda=0, mu=0, c=1, n=0)
}
\arguments{
  \item{lambda}{arrival rate}
  \item{mu}{server service rate}
  \item{c}{number of servers}
  \item{n}{number of clients in the system. Put n=0 for a standard model}
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
i_mmc <- NewInput.MMC(lambda=5, mu=10, c=2, n=0)
}
\keyword{M/M/c}

