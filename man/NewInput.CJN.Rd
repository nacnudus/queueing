% File man/NewInput.CJN.Rd
\name{NewInput.CJN}
\alias{NewInput.CJN}
\title{Define the inputs of a new one class Closed Jackson Network}
\description{
  Define the inputs of a new one class Closed Jackson Network
}

\usage{
  \method{NewInput}{CJN}(prob=NULL, n=0, z=0, ...)
}
\arguments{
  \item{prob}{It is probability transition matrix or visit ratio vector. That is, the prob[i, j] is the transition probability of node i to node j, or prob[i] is the visit ratio to node i}
  \item{n}{number of clients in the Network}
  \item{z}{think time of the client}
  \item{\dots}{a separated by comma list of nodes of i_MM1, i_MMC or i_MMInf class}
}

\details{Define the inputs of a new one class Closed Jackson Network}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{QueueingModel.i_CJN}}
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

cjn1 <- NewInput.CJN(prob, n, z, n1, n2)

}

\keyword{Closed Jackson Network}
