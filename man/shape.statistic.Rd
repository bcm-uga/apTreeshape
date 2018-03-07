\name{shape.statistic}
\alias{shape.statistic}

\title{ Computes the log of the likelihood ratio (yule/pda) }

\description{
\code{shape.statistic} computes the logarithm of the ratio of the likelihoods under the Yule model and the PDA model of the given tree. 
}

\usage{
shape.statistic(tree, norm=NULL)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
\item{norm}{ A character string equals to \code{NULL} for no normalization, \code{"yule"} for the Yule model normalization or \code{"pda"} for the pda normalization.}
}

\details{
The log of the likelihood ratio is proportional to \deqn{\sum{(\log{N(v)-1})},}{sum(log(N(v)-1)),} for all internal node \eqn{v}{v} (where \eqn{N(v)}{N(v)} is the number of internal nodes descending from the node \eqn{v}{v} ). The ratio of the likelihoods enables to build the most powerful test of the Yule model against the PDA one. (Neyman-Pearson lemma).

Under the PDA model, the log ratio has approximate Gaussian distribution of \eqn{mean \sim  2.03*n-3.545*\sqrt{n-1}}{mean ~  2.03*n-3.545*sqrt(n-1)} and \eqn{variance \sim 2.45*(n-1)*\log{n-1}}{variance ~ 2.45*(n-1)*log(n-1)}, where n is the number of tips of the tree. The Gaussian approximation is accurate for very large n (n greater than 10000(?)). The normalization of the ratio uses tabulated empirical values of variances estimated from Monte Carlo simulations. The normalization uses the formula: \deqn{variance \sim 1.570*n*\log{n}-5.674*n+3.602*\sqrt{n}+14.915}{variance ~ 1.570*n*log(n)-5.674*n+3.602*sqrt(n)+14.915}
Under the Yule model, the log ratio has approximate Gaussian distribution of \eqn{mean = 1.204*n-\log{n-1}-2}{mean = 1.204*n-log(n-1)-2} and \eqn{variance = 0.168*n-0.710}{variance = 0.168*n-0.710}, where n is the number of tips of the tree. The Gaussian approximation is accurate for n greater than 20.
}

\value{
An object of class \code{numeric} containing the shape statistic of the tree.
}

\references{
Fill, J. A. (1996), On the Distribution of Binary Search Trees under the Random Permutation Model. \emph{Random Structures and Algorithms}, \bold{8}, 1 -- 25, for more details about the normalization and proofs.
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}> \cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier François <\email{olivier.francois@imag.fr}> 
}

\examples{
data(universal.treeshape)
tree <- universal.treeshape
plot(tree)
summary(tree)

likelihood.test(tree,  model = "yule", alternative = "two.sided")
likelihood.test(tree,  model = "pda", alternative = "two.sided")

## Histogram of shape statistics for a list of Yule trees 
##      (may take some time to compute)
main="Histogram of shape statistics"; xlab="shape statistic"	
hist(sapply(rtreeshape(1000,tip.number=100,model="yule"),FUN=shape.statistic,
      norm="yule"), freq=FALSE, main=main, xlab=xlab)

## Does it fit the Gaussian distribution with mean=0 and sd=1 ?
x<-seq(-3,3,by=0.001)
lines(x,dnorm(x))
}

\keyword{ univar}

