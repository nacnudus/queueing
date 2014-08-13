% File man/VNq.o_MM1.Rd
\name{VNq.o_MM1}
\alias{VNq.o_MM1}
\title{Returns the variance of the number of customers in the queue in the M/M/1 queueing model}
\description{
  Returns the variance of the number of customers in the queue in the M/M/1 queueing model
}
\usage{
  \method{VNq}{o_MM1}(x, \dots)
}
\arguments{
  \item{x}{a object of class o_MM1}
  \item{\dots}{aditional arguments}
}

\details{Returns the variance of the number of customers in the queue in the M/M/1 queueing model}

\references{
 [Sztrik2012] Dr. Janos Sztrik (2012).\cr
\emph{Basic Queueing Theory}.\cr
University of Debrecen, Faculty of Informatics.
}

\seealso{
  \code{\link{QueueingModel.i_MM1}}.
}


\examples{
## create input parameters
i_mm1 <- NewInput.MM1(lambda=1/4, mu=1/3, n=0)

## Build the model
o_mm1 <- QueueingModel(i_mm1)

## Returns the variance
VNq(o_mm1)
}

\keyword{M/M/1}
