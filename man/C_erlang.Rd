% File man/C_erlang.Rd
\name{C_erlang}
\alias{C_erlang}
\title{Returns the probability to wait in queue bacause all servers are busy}
\description{
  Returns the probability to wait in queue bacause all servers are busy
}
\usage{
  C_erlang(c, r)
}
\arguments{
  \item{c}{numbers of servers}
  \item{r}{lambda/mu, that is, ratio of rate of arrivals and rate of service}
}

\details{Returns the probability to wait in queue bacause all servers are busy}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{B_erlang}}
}


\examples{
## two servers
C_erlang(2, 0.5/0.7)

}
\keyword{C_erlang}
