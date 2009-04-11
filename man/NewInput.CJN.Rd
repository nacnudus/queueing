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
  \item{prob}{It is probability transition matrix or visit ratio vector. That is, the prob[i, j] is the transition probability of node i to node j, or prob[i] is the visit ratio (a probability, that is, a value between 0 and 1) to node i}
  \item{n}{number of clients in the Network}
  \item{z}{think time of the client}
  \item{\dots}{a separated by comma list of nodes of i_MM1, i_MMC or i_MMInf class}
}

\details{Define the inputs of a new one class Closed Jackson Network}

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
  \code{\link{QueueingModel.i_CJN}}
}


\examples{
## See example 11.13 in reference [Sixto2004] for more details.
## create the nodes
n <- 2
n1 <- NewInput.MM1(lambda=0, mu=1/0.2, n=0)
n2 <- NewInput.MM1(lambda=0, mu=1/0.4, n=0)

## think time = 0
z <- 0

## definition of the transition probabilities
prob <- matrix(data=c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

cjn1 <- NewInput.CJN(prob, n, z, n1, n2)


## using visit ratios and service demands. See [Lazowska84] pag 117.
## E[S] cpu = 0.005, Visit cpu = 121, D cpu = E[S] cpu * Visit cpu = 0.605
cpu <- NewInput.MM1(mu=1/0.605)

## E[S] disk1 = 0.030, Visit disk1 = 70, D disk1 = E[S] disk1 * Visit disk1 = 2.1
disk1 <- NewInput.MM1(mu=1/2.1)

## E[S] disk2 = 0.027, Visit disk2 = 50, D disk2 = E[S] disk2 * Visit disk2 = 1.35
disk2 <- NewInput.MM1(mu=1/1.35)

## The visit ratios have to be set to 1. Don't put here values other than probabilities.
visit <- c(1, 1, 1)

net <- NewInput.CJN(prob=visit, n=3, z=15, cpu, disk1, disk2)

}

\keyword{Closed Jackson Network}
