% File man/WWq.o_MMCC.Rd
\name{WWq.o_MMCC}
\alias{WWq.o_MMCC}
\title{Returns the mean time spend in queue when there is queue in the M/M/c/c queueing model}
\description{
  Returns the mean time spend in queue when there is queue in the M/M/c/c queueing model
}
\usage{
  \method{WWq}{o_MMCC}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MMCC}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in queue when there is queue in the M/M/c/c queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_MMCC}}.
}


\examples{
## See example 10.12 in reference [Sixto2004] for more details.
## create input parameters
i_mmcc <- NewInput.MMCC(lambda=3, mu=0.25, c=15)

## Build the model
o_mmcc <- QueueingModel(i_mmcc)

## Returns the WWq
WWq(o_mmcc)
}

\keyword{M/M/c/c}
