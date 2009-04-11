% File man/Lq.Rd
\name{Lq}
\alias{Lq}
\title{Returns the mean number of clients in the queue in a queueing model}
\description{
  Returns the mean number of clients in the queue in a queueing model
}
\usage{
  Lq(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of clients in the queue in a queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{Lq.o_MM1}}\cr
  \code{\link{Lq.o_MMC}}\cr
  \code{\link{Lq.o_MM1K}}\cr
  \code{\link{Lq.o_MMCK}}\cr
  \code{\link{Lq.o_MM1KK}}\cr
  \code{\link{Lq.o_MMCKK}}\cr
  \code{\link{Lq.o_MMCC}}\cr
  \code{\link{Lq.o_MMCKM}}\cr
  \code{\link{Lq.o_MMInf}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the Lq
Lq(o_mm1)
}

\keyword{M/M/1}
