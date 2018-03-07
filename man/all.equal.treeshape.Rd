\name{all.equal.treeshape}
\alias{all.equal.treeshape}

\title{ Compare two objects of class treeshape }

\description{
This function makes a global comparison of two phylogenetic trees.
}

\usage{
\method{all.equal}{treeshape}(target, current, names=FALSE, height=FALSE, \dots)
}

\arguments{
\item{target}{An object of class \code{"treeshape"}.}
\item{current}{An object of class \code{"treeshape"}.}
\item{names}{An object of class \code{"logical"}, checking if the names of the tips should be tested. If \code{FALSE} (default), the names of the tips are not compared.}
\item{height}{An object of class \code{"logical"}, checking if the heights of the nodes should be tested. If \code{FALSE} (default), the height of internal nodes are not compared.}
\item{...}{further arguments passed to or from other methods.}
}

\details{
This function is meant to be an adaptation of the generic function \code{all.equal} for the comparison of phylogenetic trees. A phylogenetic tree can have many different representations. Permutations between the left and the right daughter clade of a node do not change the corresponding phylogeny, and \code{all.equal.treeshape} returns \code{TRUE} on two permutated trees.
}

\value{
Returns the logical \code{TRUE} if the tree objects are similar up to a permutation of their tips. Otherwise, it returns \code{FALSE}. Heights and labels can be taken into account.
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}>\cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{
\code{\link[base]{all.equal}} for the generic R function.
}

\examples{
  
## Trees with permutations 
data(carnivora.treeshape)
tree=carnivora.treeshape
tree$merge[8,]=c(tree$merge[8,2],tree$merge[8,1])
all.equal(tree, carnivora.treeshape)
 
## Trees with different heights
merge=matrix(NA, 3, 2)
merge[,1]=c(-3,-1,2); merge[,2]=c(-4,-2,1);tree1=treeshape(merge)
merge[,1]=c(-1,-3,1); merge[,2]=c(-2,-4,2);tree2=treeshape(merge)
  
plot(tree1, tree2)
all.equal(tree1, tree2)
all.equal(tree1, tree2, height=TRUE)

## Trees with different names
tree3=treeshape(tree1$merge, c("a", "b", "c", "d"))
tree4=treeshape(tree1$merge, c("1", "2", "3", "4"))
plot(tree3, tree4)
all.equal(tree3, tree4)
all.equal(tree3, tree4, names=TRUE)

}

\keyword{ manip }
