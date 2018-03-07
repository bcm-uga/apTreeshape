\name{summary.treeshape}
\alias{summary.treeshape}

\title{ Print a summary of an object of class "treeshape"}

\usage{ 
\method{summary}{treeshape}(object, \dots)
}

\arguments{
\item{object}{an object of class \code{"treeshape"}.}
\item{...}{further arguments passed to or from other methods.}
}

\description{
This function prints a compact summary of a phylogenetic tree of class \code{"treeshape"})
}
\details{
\code{summary.treeshape} prints the following information: the number of tips of the tree, its Colless' index, and the expected values and standard deviations of the Colless' index under the PDA and Yule models.The expected value of the Colless' index under the Yule model is given according to the formula: \eqn{n*\log{n}+n*(\gamma-1-\log{2})}{n*log(n)+n*(gamma-1-log(2))} where n is the number of tips of the tree and \eqn{\gamma}{gamma} the Euler's constant. The standard deviation under the Yule model is given by: \eqn{\sqrt{(3-\frac{\pi^2}{6}-\log{2})*n}}{sqrt((3-pi^2/6-log(2))*n)}. The expected value of the Colless' index under the PDA model is given according to the formula: \eqn{\sqrt{\pi}*n^{3/2}}{sqrt(pi)*n^(3/2)}. The standard deviation under the PDA model is given by: \eqn{\sqrt{\frac{10}{3}-\pi}*n^{3/2}}{sqrt(10/3-pi)*n^(3/2)}.
}

\value{
A NULL value is returned, the results are simply printed.
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}> \cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier François <\email{olivier.francois@imag.fr}>
}

\seealso{ 
\code{\link[base]{summary}}\cr
\code{\link{colless}} for more informations about the expected values under the Yule and PDA models.\cr
}

\examples{

## Summary of a PDA tree with 100 tips.
summary(rpda(100)) 
## Note that the standard deviation is very large. 

## Summary of a Yule tree with 100 tips.
summary(ryule(100)) 
## The standard deviation under the Yule model is much smaller than under 
##      the PDA model.

## Summary of the HIV tree.
data(hivtree.treeshape)
summary(hivtree.treeshape) 
## The HIV tree is much closer from the Yule model than from the PDA model. 
}

\keyword{manip}
