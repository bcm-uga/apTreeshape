\name{spectrum.treeshape}
\alias{spectrum.treeshape}

\title{ Compute the spectrum of a tree}

\description{
This function returns a sequence containing the number of subtrees of size n, n-1, ..., 3, 2 where n is the size of the tree. The 'k'th element of the sequence is the number of subtrees of size n-k+1 in the tree, where n is the number of tips of the tree. 
}

\usage{
spectrum.treeshape(tree)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
}

\value{
A sequence of size n-1 (where n is the number of tips of the tree) containing the number of subtrees of each size.
\code{spectrum.treeshape(tree)[1]} is the number of subtrees with n tips (equal to 1).
\code{spectrum.treeshape(tree)[n-1]} is the number of cherries of the tree (subtrees with 2 tips).
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{
\code{\link{smaller.clade.spectrum}} 
}

\examples{

## A random Yule tree with 30 tips
tr<-rtreeshape(n=1,tip.number=30,model="yule")
tr<-tr[[1]]
spectre=spectrum.treeshape(tr)
spectre
  
## Number of cherries of the tree : nrow(tr$merge)==29
spectre[29]
}

\keyword{univar}
