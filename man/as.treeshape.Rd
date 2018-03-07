\name{as.treeshape}
\alias{as.treeshape}
\alias{as.treeshape.phylo}

\title{ Conversion among tree objects }

\description{
\code{as.treeshape} is a generic function which converts an tree object into an object of class \code{"treeshape"}. There is currently a single method for this generic function. 
}
  
\usage{
\method{as.treeshape}{phylo}(x, model, p, ...)
}

\details{
\code{as.treeshape} can convert trees that are not binary. When trying to convert a tree with polytomies, this function may either reject the tree (if \code{model=NULL}) or simulate the tree. The polytomy is replaced by a randomized subtree with n tips where n is the size of the polytomy. The subtree is simulated using the PDA, Yule, Aldous or biased model.
}

\arguments{
\item{x}{An object to be converted into another class. Currently, it must be an object of class \code{"phylo"}.}
\item{model}{The model to use when the tree to convert is not binary. If \code{NULL} (default), the tree is not converted. One of \code{"biased"}, \code{"pda"}, \code{"aldous"} or \code{"yule"} character string.}
\item{p}{The parameter for the model \code{"biased"}.}
\item{...}{Further arguments to be passed to or from other methods.}
}

\value{
An object of class \code{"treeshape"} or an object of classes \code{"treeshape"} and \code{"randomized.treeshape"} if the original tree was not binary.
}

\author{
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{ 
\code{\link{as.phylo.treeshape}}\cr
}

\examples{

library(ape)
data(bird.orders)
## Data set from APE
plot(bird.orders)
  
## "treeshape" conversion
tree=as.treeshape(bird.orders)
plot(tree)
summary(tree)
 
}

\keyword{manip}

