% File man/Wc.Rd
\name{Wc}
\alias{Wc}
\title{Returns the vector with each class mean time spend on a multiclass queueing network}
\description{
  Returns the vector with each class mean time spend on a multiclass queueing network
}
\usage{
  Wc(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCON, o_MCCN}
  \item{\dots}{aditional arguments}
}

\details{Returns the vector with each class mean time spend on a multiclass queueing network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.\cr
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{Wc.o_MCON}}\cr
  \code{\link{Wc.o_MCCN}}
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

classes <- 2
vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vService <- matrix(data=c(1/10, 1/3, 2/5, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(classes, vLambda, nodes, vType, vVisit, vService)

# Build the model
o_mcon1 <- QueueingModel(i_mcon1)

Wc(o_mcon1)

}
\keyword{MultiClass Queueing Network}
