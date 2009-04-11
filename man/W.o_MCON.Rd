% File man/W.o_MCON.Rd
\name{W.o_MCON}
\alias{W.o_MCON}
\title{Returns the mean time spend in a MultiClass Open Network}
\description{
  Returns the mean time spend in a MultiClass Open Network
}
\usage{
  \method{W}{o_MCON}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MCON}
  \item{\dots}{aditional arguments}
}

\details{Returns the mean time spend in a MultiClass Open Network}

\references{
  [Lazowska84] Edward D. Lazowska, John Zahorjan, G. Scott Graham, and Kenneth C. Sevcik (1984).\cr
\emph{Quantitative System Performance: Computer System Analysis Using Queueing Network Models}.\cr
Prentice-Hall, Inc., Englewood Cliffs, New Jersey
}

\seealso{
  \code{\link{QueueingModel.i_MCON}}.
}


\examples{
## See example in pag 138 in reference [Lazowska84] for more details.

vLambda <- c(3/19, 2/19)
nodes <- 2
vType <- c("Q", "Q")
vVisit <- matrix(data=c(10, 9, 5, 4), nrow=2, ncol=2, byrow=TRUE)
vMu <- matrix(data=c(10, 3, 5/2, 1), nrow=2, ncol=2, byrow=TRUE)

i_mcon1 <- NewInput.MCON(vLambda, nodes, vType, vVisit, vMu)

# Build the model
o_mcon1 <- QueueingModel(i_mcon1)

W(o_mcon1)

}

\keyword{MultiClass Open Network}
