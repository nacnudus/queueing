% File man/CheckInput.i_CJN.Rd
\name{CheckInput.i_CJN}
\alias{CheckInput.i_CJN}
\title{Check the input params of a Closed Jackson Network}
\description{
  Check the input params of a Closed Jackson Network
}
\usage{
  \method{CheckInput}{i_CJN}(x, \dots)
}
\arguments{
  \item{x}{a object of class i_CJN}
  \item{\dots}{aditional arguments}
}

\details{Check the input params of a Closed Jackson Network. The inputs params are created calling previously the \link{NewInput.CJN}}

\references{
  [Sixto2004] Sixto Rios Insua, Alfonso Mateos Caballero, M Concepcion Bielza Lozoya, Antonio Jimenez Martin (2004).\cr
\emph{Investigacion Operativa. Modelos deterministicos y estocasticos}.\cr
Editorial Centro de Estudios Ramon Areces.
}

\seealso{
  \code{\link{NewInput.CJN}}
}


\examples{
## See example 11.13 in reference [Sixto2004] for more details.
## create the nodes
n <- 2
n1 <- NewInput.MM1(lambda=0, mu=1/0.2, n=0)
n2 <- NewInput.MM1(lambda=0, mu=1/0.4, n=0)

# think time = 0
z <- 0

# operational value
operational <- FALSE

# definition of the transition probabilities
prob <- matrix(data=c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

cjn1 <- NewInput.CJN(prob, n, z, operational, n1, n2)

CheckInput(cjn1)

}
\keyword{Closed Jackson Network}
