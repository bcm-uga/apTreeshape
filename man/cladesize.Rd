\name{cladesize}
\alias{cladesize}

\title{ Compute the number of children of a randomly chosen node }

\description{
\code{cladesize} takes a random internal node in a tree, and computes its number of descendants (clade size).
}
\usage{
cladesize(tree)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}. }
}

\details{ 
This function can be used to check whether a tree fits the Yule or the PDA models. Under the Yule model, the probability distribution of the random clade size is equal to \deqn{P(Kn=k)=\frac{2n}{(n-1)k(k+1)}}{P(Kn=k)=2*n/((n-1)*k(k+1))} for \eqn{k = 2, 3, \ldots, n-1} and \deqn{P(Kn=n)=\frac{1}{n-1}}{P(Kn=n)=1/(n-1)} (where \eqn{n}{n} is the number of tips of the tree and \eqn{Kn}{Kn} is the number of descendents of an internal node of the tree). Under the PDA model, the asymptotic distribution (when the number of tips grows to infinity) of the random clade size is equal to \deqn{P(K=k+1)=\frac{{2k\choose k}}{(k+1)(2^k)^2}}{P(K=k+1)=choose(2*k,k)/((k+1)*(2^k)^2)}.
}

\value{
An object of class \code{numeric} representing the clade size of a random node of a tree.
}

\references{
Blum, M., Francois, O. and Janson, S. The mean, variance and limiting distribution of two statistics sensitive to phylogenetic tree balance; manuscript available from \cr  \url{http://www-timc.imag.fr/Olivier.Francois/bfj.pdf}.
}

\author{
Michael Blum <\email{michael.blum@imag.fr}> \cr 
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Oliver François <\email{olivier.francois@imag.fr}>
}
 
\examples{

# Histogram of random clade sizes 
main="Random clade sizes for random generated trees"
xlabel="clade size"
hist(sapply(rtreeshape(100,tip.number=40,model="yule"),FUN=cladesize),
      freq=FALSE,main=main,xlab=xlabel)
}

\keyword{univar}