% File man/Throughput.o_MMCKK.Rd
\name{Throughput.o_MMCKK}
\alias{Throughput.o_MMCKK}
\title{Throughput of a M/M/c/K/K queueing model}
\description{
  Returns the throughput of a M/M/c/K/K queueing model
}
\usage{
  \method{Throughput}{o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a M/M/c/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.MMCKK}}, \code{\link{CheckInput.i_MMCKK}}, \code{\link{QueueingModel.i_MMCKK}}
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12)

## build the model
o_mmckk <- QueueingModel(i_mmckk)

## Throughput
Throughput(o_mmckk)

}
\keyword{M/M/c/K/K}
