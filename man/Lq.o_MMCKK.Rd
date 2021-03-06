% File man/Lq.o_MMCKK.Rd
\name{Lq.o_MMCKK}
\alias{Lq.o_MMCKK}
\title{Returns the mean number of customers in the queue in the M/M/c/K/K queueing model}
\description{
  Returns the mean number of customers in the queue in the M/M/c/K/K queueing model
}
\usage{
  \method{Lq}{o_MMCKK}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCKK}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of customers in the queue in the M/M/c/K/K queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCKK}}.
}


\examples{
## create input parameters
i_mmckk <- NewInput.MMCKK(lambda=8, mu=2, c=5, k=12, method=0)

## Build the model
o_mmckk <- QueueingModel(i_mmckk)

## Returns the Lq
Lq(o_mmckk)
}

\keyword{M/M/c/K/K}
