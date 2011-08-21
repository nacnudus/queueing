% File man/ROk.Rd
\name{ROk}
\alias{ROk}
\title{Reports a vector with each node (server) use of a queueing network}
\description{
  Reports a vector with each node (server) use of a queueing network
}
\usage{
 ROk(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCON, o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{Reports a vector with each node (server) use of a queueing network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos CaballeROk, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial CentROk de Estudios Ramon Areces.\cr\cr

[Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{ROk.o_OJN}}\cr
  \code{\link{ROk.o_CJN}}\cr
  \code{\link{ROk.o_MCON}}\cr
  \code{\link{ROk.o_MCCN}}
}


\examples{
## See example in pag 142 in reference [Lazowska84] for more details.

classes <- 2
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_MCCN1 <- NewInput.MCCN(classes, vNumber, vThink, nodes, vType, vVisit, vMu)

# Build the model
o_MCCN1 <- QueueingModel(i_MCCN1)

ROk(o_MCCN1)

}
\keyword{Queueing Networks}
