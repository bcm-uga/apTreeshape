\name{index.test}
\alias{colless.test}
\alias{sackin.test}

\title{ Perform a test on the Yule or PDA hypothesis based on the Colless or the Sackin statistic}

\description{
\code{colless.test} performs a test based on the Colless' index on tree data for the Yule or PDA model hypothesis. \cr
\code{sackin.test} does the same with the Sackin's index. \cr
}

\usage{
colless.test(tree, model = "yule", alternative = "less", n.mc = 500)
sackin.test(tree, model = "yule", alternative = "less", n.mc = 500)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.  }
\item{model}{ The null hypothesis of the test. One of  \code{"yule"} (default) or \code{"pda"}. }
\item{alternative}{ A character string specifying the alternative hypothesis of the test. One of \code{"less"} (default) or \code{"greater"}. }
\item{n.mc}{ An integer representing the number of random trees to be generated and required to compute a p-value from a Monte Carlo method. }
}

\details{
A test on tree data that either rejects the Yule or the PDA models.
This test is based on a Monte Carlo estimate of the p-value. Replicates are generated under the Yule or PDA models, and their Colless' (Sackin's) indices are computed. The empirical distribution function of these statistics is then computed thanks to the "ecdf" R function. The p-value is then deduced from its quantiles. The less balanced the tree is and the larger its Colless's (Sackin's) index. The alternative "less" should be used to test whether the tree is more balanced (less unbalanced) than predicted by the null model. The alternative "greater" should be used to test whether the tree is more unbalanced than predicted by the null model. The computation of p-values may take some  time depending on the number of replicates |n.mc| and the size of the simulated tree.
}

\value{
\item{model}{the null model.}
\item{statistic}{the test statistic.}
\item{p.value}{the p-value of the test.}
\item{alternative}{the alternative hypothesis of the test.}
}

\references{ 
Mooers, A. O., Heard, S. B. (Mar., 1997) Inferring Evolutionnary Process from Phylogenetic Tree Shape. \emph{The Quarterly Review of Biology}, \bold{72}, 31 -- 54.
  
Blum, M., Francois, O. and Janson, S. The mean, variance and limiting distribution of two statistics sensitive to phylogenetic tree balance; manuscript available from \url{http://www-timc.imag.fr/Olivier.Francois/bfj.pdf}.
}

\author{
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr }>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{ 
\code{\link{subtree.test}}\cr
\code{\link{colless}}\cr
\code{\link{sackin}}\cr
}

\examples{

## Test on a randomly generated Yule tree with 30 tips
a<-rtreeshape(1,30,model="yule")
a<-a[[1]]
  
## Is it more balanced than a Yule tree ?
colless.test(a,alternative="less",model="yule")
## Is it less balanced than a PDA tree ?
colless.test(a,model="pda",alternative="greater")
 
## Test on the phylogenetic tree hiv.treeshape: is it more balanced than 
##      predicted by the Yule model?
data(hivtree.treeshape)
## The tree looks compatible with the null hypothesis
colless.test(hivtree.treeshape, alternative="greater", model="yule")
 
## What happen when we look at the top the tree?
colless.test(cutreeshape(hivtree.treeshape, 160, "top"),
      alternative="greater", model="yule")
colless.test(cutreeshape(hivtree.treeshape, 160, "top"), 
      alternative="greater", model="pda")

## Test with the Sackin's index: is the HIV tree less balanced than 
##      predicted by the PDA model?
sackin.test(hivtree.treeshape,alternative="greater",model="pda") 
## The p.value equals to 1...
}
  
\keyword{htest}
