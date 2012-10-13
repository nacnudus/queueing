% File man/Throughputk.o_MCCN.Rd
\name{Throughputk.o_MCCN}
\alias{Throughputk.o_MCCN}
\title{Reports a vector with each node (server) throughput of a MultiClass Closed Network}
\description{
  Reports a vector with each node (server) throughput of a MultiClass Closed Network
}
\usage{
  \method{Throughputk}{o_MCCN}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{
Reports a vector with each node (server) throughput of a MultiClass Closed Network
}

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
  \code{\link{NewInput.MCCN}}, \code{\link{CheckInput.i_MCCN}}, \code{\link{QueueingModel.i_MCCN}}
}


\examples{
## See example in pag 142 in reference [Lazowska84] for more details.

classes <- 2
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vService <- matrix(data=c(1/10, 1/3, 2/5, 1), nrow=2, ncol=2, byrow=TRUE)

i_MCCN1 <- NewInput.MCCN(classes, vNumber, vThink, nodes, vType, vVisit, vService)

# Build the model
o_MCCN1 <- QueueingModel(i_MCCN1)

Throughputk(o_MCCN1)

}
\keyword{MultiClass Closed Network}

