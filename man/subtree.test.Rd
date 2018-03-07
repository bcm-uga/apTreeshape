\name{subtree.test}
\alias{subtree.test}

\title{ Test the Yule or PDA hypothesis  }

\description{
\code{subtree.test} tests the likelihood of the Yule or the PDA hypothesis and computes the p-value of the test. The test is based on the number of subtrees of a given size in the tree.
}

\usage{
subtree.test(tree, size , alternative = "two.sided")
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
\item{size}{ The size of the subtrees to test for.}
\item{alternative}{ The alternative hypothesis of the test. It can be \code{"two.sided"} (default), \code{"less"} or \code{"greater"}.}
}

\details{
See references for the mathematical details of the test. It uses a Gaussian approximation to compute the p-value.
}

\value{
A list containing the following arguments :
\item{statistic}{the value of the statistic used in the test}
\item{p.value }{the p-value of the test}
\item{alternative}{the alternative hypothesis used for the test}
}

\references{
Blum, M. and François, O. (2005) Minimal clade size and external branch length under the neutral coalescent \emph{Advances in Applied Probability} (2005), \bold{37.3}. 
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}


\seealso{ 
\code{\link{sackin.test}} \cr
\code{\link{colless.test}} 
}

\examples{

## Generate a random pda tree with 50 tips
tr<-rtreeshape(n=1,tip.number=50,model="pda")
tr<-tr[[1]]

## Test the yule hypothesis, using subtrees of size 2 (Cherries), 
##      with the alternative hypothesis "less"
subtree.test(tr,size=2,alternative="less")
}

\keyword{htest}
