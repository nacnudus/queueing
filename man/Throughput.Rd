% File man/Throughput.Rd
\name{Throughput}
\alias{Throughput}
\title{Throughput of a queueing model (or network)}
\description{
  Returns the throughput of a queueing model (or network)
}
\usage{
 Throughput(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1, o_MMC, o_MM1K, o_MMCK, o_MM1KK, o_MMCKK, o_MMCC, o_MMCKM, o_MMInfKK, o_MMInf, o_OJN, o_CJN, o_MCON, o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{Returns the throughput of a queueing model (or network)}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
\cr
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{Throughput.o_MM1}}\cr
  \code{\link{Throughput.o_MMC}}\cr
  \code{\link{Throughput.o_MM1K}}\cr
  \code{\link{Throughput.o_MMCK}}\cr
  \code{\link{Throughput.o_MM1KK}}\cr
  \code{\link{Throughput.o_MMCKK}}\cr
  \code{\link{Throughput.o_MMCC}}\cr
  \code{\link{Throughput.o_MMCKM}}\cr
  \code{\link{Throughput.o_MMInfKK}}\cr
  \code{\link{Throughput.o_MMInf}}\cr
  \code{\link{Throughput.o_OJN}}\cr
  \code{\link{Throughput.o_CJN}}\cr
  \code{\link{Throughput.o_MCON}}\cr
  \code{\link{Throughput.o_MCCN}}
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Throughput
Throughput(o_mm1)

}
\keyword{Queueing Models}
