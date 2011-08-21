% File man/Throughputc.o_MCON.Rd
\name{Throughputc.o_MCON}
\alias{Throughputc.o_MCON}
\title{Reports a vector with each class throughput in a MultiClass Open Network}
\description{
  Reports a vector with each class throughput in a MultiClass Open Network
}
\usage{
  \method{Throughputc}{o_MCON}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCON}
  \item{\dots}{aditional arguments}
}

\details{
  Reports a vector with each class throughput in a MultiClass Open Network
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
  \code{\link{NewInput.MCON}}, \code{\link{CheckInput.i_MCON}}, \code{\link{QueueingModel.i_MCON}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

classes <- 2
vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(classes, vLambda, nodes, vType, vVisit, vMu)

# Build the model
o_mcon1 <- QueueingModel(i_mcon1)

Throughputc(o_mcon1)

}
\keyword{MultiClass Open Network}

