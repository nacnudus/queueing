% File man/ROck.o_MCMN.Rd
\name{ROck.o_MCMN}
\alias{ROck.o_MCMN}
\title{Reports a matrix with the use of class i in each node (server) j in a MultiClass Mixed Network}
\description{
  Reports a matrix with the use of class i in each node (server) j in a MultiClass Mixed Network
}
\usage{
  \method{ROck}{o_MCMN}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCMN}
  \item{\dots}{aditional arguments}
}

\details{Reports a matrix with the use of class i in each node (server) j in a }

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{QueueingModel.i_MCMN}}.
}


\examples{
## See example in pag 147 in reference [Lazowska84] for more details.

classes <- 4
vLambda <- c(1, 1/2)
vNumber <- c(1, 1)
vThink <- c(0, 0)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=1, nrow=4, ncol=2)
vService <- matrix(data=c(1/4, 1/2, 1/2, 1, 1/6, 1, 1, 4/3), nrow=4, ncol=2)

i_mcmn1 <- NewInput.MCMN(classes, vLambda, vNumber, vThink, nodes, vType, vVisit, vService)

# Build the model
o_mcmn1 <- QueueingModel(i_mcmn1)

ROck(o_mcmn1)

}
\keyword{MultiClass Mixed Network}
