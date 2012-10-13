% File man/Engset.Rd
\name{Engset}
\alias{Engset}
\title{Returns the probability that all servers are busy}
\description{
  Returns the probability that all servers are busy
}
\usage{
  Engset(k=1, c=0, r=0)
}
\arguments{
  \item{k}{numbers of users}
  \item{c}{numbers of servers}
  \item{r}{lambda/mu, that is, ratio of rate of arrivals and rate of service}
}

\details{Returns the probability of blocking in a finite source model}

\seealso{
  \code{\link{B_erlang}}
}


\examples{
## three users, two servers
Engset(3, 2, 0.5/0.7)
}
\keyword{Engset}
