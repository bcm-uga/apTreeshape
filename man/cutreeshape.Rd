\name{cutreeshape}
\alias{cutreeshape}

\title{ Cut objects of class "treeshape" }

\description{
Prunes or cuts an object of class \code{"treeshape"} from a specificed internal node, either by specifying a top or bottom direction. This function returns either the top part or the bottom part of a tree. 
}

\usage{
cutreeshape(tree, node, type)
}

\arguments{
\item{tree}{An object of class \code{"treeshape"}.}
\item{node}{An integer representing the node at which the tree will be pruned or cut. \code{node} should be in the range \eqn{1, \ldots, treesize-1}{1:(treesize-1)}.}
\item{type}{A character string equals to either \code{"top"} or \code{"bottom"}.}
}

\details{

If the \code{type} argument is "top", the tree is pruned from \code{node}. The resulting tips correspond to the ancestral branches present at the same height as the given node. New tip labels are assigned to the tips.\cr
  
If the \code{type} specified is "bottom", the subtree under \code{node} is returned. The tips are not renamed (they keep their former names) and the specified \code{node} is the root of the new tree.
  
}
\value{
An object of class \code{"treeshape"}
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag.fr}>
}

\seealso{
\code{\link{tipsubtree}}
}

\examples{

## Data set provided with the library. Type help(cytochromc) for more infos.
data(carnivora.treeshape)  
data(hivtree.treeshape)

## Examples of "bottom" cutting:
bottom.tree=cutreeshape(carnivora.treeshape, 3, "bottom")
plot(carnivora.treeshape, bottom.tree)
bottom.tree=cutreeshape(carnivora.treeshape, 8, "bottom")
plot(carnivora.treeshape, bottom.tree)
  
## Examples of "top" pruning:
top.tree=cutreeshape(hivtree.treeshape, 158, "top")
plot(hivtree.treeshape, top.tree)
}

\keyword{manip}
