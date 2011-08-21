% File man/NewInput.CJN.Rd
\name{NewInput.CJN}
\alias{NewInput.CJN}
\title{Define the inputs of a new one class Closed Jackson Network}
\description{
  Define the inputs of a new one class Closed Jackson Network
}

\usage{
  \method{NewInput}{CJN}(prob=NULL, n=0, z=0, operational=FALSE, ...)
}
\arguments{
  \item{prob}{It is probability transition matrix or visit ratio vector. That is, the prob[i, j] is the transition probability of node i to node j, or prob[i] is the visit ratio (a probability, that is, a value between 0 and 1) to node i. Also, the visit ratio can express the number of times that a client visits the queueing center, in a more operational point of view. See the parameter operational}
  \item{n}{number of clients in the Network}
  \item{z}{think time of the client}
  \item{operational}{If prob is a vector with the visit ratios, operational equal to FALSE gives to the visit ratio a probability meaning, that is, as the stacionary values of the imbedded markov chain. If operational is equal to TRUE, the operational point of view is used: it is the number of visits that the same client makes to a node.}
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

## operational value
operational <- FALSE

## definition of the transition probabilities
prob <- matrix(data=c(0.5, 0.5, 0.5, 0.5), nrow=2, ncol=2, byrow=TRUE)

cjn1 <- NewInput.CJN(prob, n, z, operational, n1, n2)


## using visit ratios and service demands. See [Lazowska84] pag 117.
## E[S] cpu = 0.005, Visit cpu = 121, D cpu = E[S] cpu * Visit cpu = 0.605
cpu <- NewInput.MM1(mu=1/0.005)

## E[S] disk1 = 0.030, Visit disk1 = 70, D disk1 = E[S] disk1 * Visit disk1 = 2.1
disk1 <- NewInput.MM1(mu=1/0.030)

## E[S] disk2 = 0.027, Visit disk2 = 50, D disk2 = E[S] disk2 * Visit disk2 = 1.35
disk2 <- NewInput.MM1(mu=1/0.027)

## The visit ratios.
visit <- c(121, 70, 50)

operational <- TRUE

net <- NewInput.CJN(prob=visit, n=3, z=15, operational, cpu, disk1, disk2)

}

\keyword{Closed Jackson Network}
