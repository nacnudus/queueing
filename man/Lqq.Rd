% File man/Lqq.Rd
\name{Lqq}
\alias{Lqq}
\title{Returns the mean number of clients in queue when there is queue in a queueing model}
\description{
  Returns the mean number of clients in queue when there is queue in a queueing model
}
\usage{
  Lqq(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean number of clients in queue when there is queue in a queueing model}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{Lqq.o_MM1}}\cr
  \code{\link{Lqq.o_MMC}}\cr
  \code{\link{Lqq.o_MM1K}}\cr
  \code{\link{Lqq.o_MMCK}}\cr
  \code{\link{Lqq.o_MM1KK}}\cr
  \code{\link{Lqq.o_MMCKK}}\cr
  \code{\link{Lqq.o_MMCC}}\cr
  \code{\link{Lqq.o_MMCKM}}\cr
  \code{\link{Lqq.o_MMInfKK}}\cr
  \code{\link{Lqq.o_MMInf}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the Lqq
Lqq(o_mm1)
}

\keyword{M/M/1}
