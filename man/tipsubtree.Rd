\name{tipsubtree}
\alias{tipsubtree}

\title{Extract a subtree that contains pre-specified tip names or labels}

\description{
\code{tipsubtree} returns an object of class \code{"treeshape"} that contains pre-specified tip names or labels. The name of the tips are conserved. It extracts the smallest tree that contains the common ancestors of the given tips.
}

\usage{
tipsubtree(tree, tips, numeric=FALSE)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
\item{tips}{ A vector that contains the names of the tips one want to keep. Warning, the names are case-sensitive.}
\item{numeric}{ An object of class \code{"logical"}. If \code{FALSE}, the vector \code{tips} contains the names of the tips. If \code{TRUE}, it contains the labels of the tips.}
}

\value{
An object of class \code{"treeshape"}
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}> \cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier François <\email{olivier.francois@imag.fr}>
}

\seealso{
\code{\link{cutreeshape}}
}

\examples{
  
## The universal tree of life provided in the data sets.
data(universal.treeshape)
  
## One might want to extract the tree containing the Animals, the Plants,
##      the Aquifex and the Microsporidia
tree1<-tipsubtree(universal.treeshape,tips=c("Animals", "Aquifex", 
      "Microsporidia", "Plants"))
plot(universal.treeshape, tree1)

## Labels that do not appear in the tree are ignored
tree2<-tipsubtree(universal.treeshape,tips=c("Human", "Animals", "Aquifex", 
      "Microsporidia", "Plants"))
plot(universal.treeshape, tree2)
  
tree3<-tipsubtree(universal.treeshape,tips=c(1,3,7), numeric=TRUE)
plot(universal.treeshape, tree3)
}

\keyword{manip}
