% File man/summary.o_CJN.Rd
\name{summary.o_CJN}
\alias{summary.o_CJN}
\title{Reports the results of a Closed Jackson Network}
\description{
  Reports the results of a Closed Jackson Network
}
\usage{
  \method{summary}{o_CJN}(object, \dots)
}
\arguments{
  \item{object}{a object of class o_CJN}
  \item{\dots}{aditional arguments}
}

\details{generates a report of the queueing network received as parameter}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_CJN}}.
}


\examples{
## See example 11.13 in reference [Sixto2004] for more details.
## create the nodes
n <- 2
n1 <- NewInput.MM1(lambda=0, mu=1/0.2, n=0)
n2 <- NewInput.MM1(lambda=0, mu=1/0.4, n=0)

# think time = 0
z <- 0

# definition of the transition probabilities
prob <- matrix(data=c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

# Define a new input
cjn1 <- NewInput.CJN(prob, n, z, n1, n2)

# Check the inputs and build the model
m_cjn1 <- QueueingModel(cjn1)

summary(m_cjn1)

}
\keyword{Closed Jackson Network}
