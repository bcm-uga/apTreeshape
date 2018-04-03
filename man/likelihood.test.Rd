\name{likelihood.test}
\alias{likelihood.test}

\title{ Test the Yule model vs PDA (uniform) model.}

\description{
\code{likelihood.test} uses the function \code{shape.statistic} to test the Yule model against the PDA model. The test is based on a Gaussian approximation for the log-ratio of likelihoods.
}

\usage{
likelihood.test(tree, model = "yule", alternative="two.sided")
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"} on which the test is performed. }
\item{model}{ The null hypothesis of the test. It must be equal to one of the two character strings \code{"yule"} or \code{"pda"}. }
\item{alternative}{ A character string specifying the alternative hypothesis of the test. Must be one of \code{"two.sided"} (default), \code{"less"} or \code{"greater"}. }
}

\details{ 
A test on tree data that either rejects the Yule of the PDA model. The test is based on the ratio of the likelihood of the PDA model to the likelihood of the Yule model (shape.statistic). The less balanced the tree is the larger its shape statistic is. The alternative "less" should be used to test whether the tree is less unbalanced than predicted by the null model. The alternative "greater" should be used to test whether the tree is more unbalanced than predicted by the null model. \cr

Under the Yule model, the test statistic has approximate Gaussian distribution of \eqn{mean = 1.204*n-\log{n-1}-2}{mean = 1.204*n-log(n-1)-2} and \eqn{variance = 0.168*n-0.710}{variance = 0.168*n-0.710}, where \eqn{n}{n} is the number of tips of the tree. The Gaussian approximation is accurate for n greater than 20. \cr

Under the PDA model, the test statistic has approximate Gaussian distribution of \eqn{mean \sim  2.03*n-3.545*\sqrt{n-1}}{mean ~  2.03*n-3.545*sqrt(n-1)} and \eqn{variance \sim 2.45*(n-1)*\log{n-1}}{variance ~ 2.45*(n-1)*log(n-1)}, where \eqn{n}{n} is the number of tips of the tree. The Gaussian approximation is however accurate for very large n (n greater than 10000(?)). The values of the means and variances have been obtained from an analogy with binary search tree models in computer science. \cr

The function includes corrections for small sizes under the PDA model, and uses empirical values of variances estimated throught Monte Carlo replicates as follows \deqn{variance \sim 1.570*n*\log{n}-5.674*n+3.602*\sqrt{n}+14.915}{variance ~ 1.570*n*log(n)-5.674*n+3.602*sqrt(n)+14.915} \cr

}

\value{
\code{likelihood.test} returns a list which includes:
\item{model}{the null model used by the test}
\item{statistic}{the test statistic}
\item{p.value}{the p.value of the test}
\item{alternative}{the alternative hypothesis of the test}
}

\references{ 
Fill, J. A. (1996), On the Distribution of Binary Search Trees under the Random Permutation Model. \emph{Random Structures and Algorithms}, \bold{8}, 1 -- 25.\cr
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}> \cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{ \code{\link{shape.statistic}} }

\examples{

## Generate a Yule tree with 150 tips. Is it likely to be fitting the PDA model?
likelihood.test(ryule(150),model="pda") 
## The p.value is close from 0. We reject the PDA hypothesis.

## Test on the Carnivora tree: is it likely to be fitting the Yule model?
data(carnivora.treeshape)
likelihood.test(carnivora.treeshape) 
## The p.value is high, so it's impossible to reject the Yule hypothesis.
}

\keyword{ htest }
