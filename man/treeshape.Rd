\name{treeshape}
\alias{treeshape}

\title{ Builds an object of class treeshape }

\description{
\code{treeshape} builds a tree of class \code{"treeshape"} from a n*2 matrix, where n is the number of internal nodes of the tree. There are no informations about the heights of the branches in an object of class \code{"treeshape"}. Formally, a "tree shape" is a phylogenetic tree where the label of the tips are ignored. Here, the label of the tips can be kept or ignored. If a \code{names} vector is provided, then the names of species are attached to the tips of the tree. Otherwise, tips are simply labeled with their numbers in the tree. Trees of class \code{"treeshape"} are always binary.   
}

\usage{
treeshape(nodes, names)
}

\arguments{
\item{nodes}{ \code{nodes} is a n*2 matrix containing the node structure of the tree.}
\item{names}{ \code{names} is a vector which contains the names of the tips.}
}

\details{
A tree of class \code{"treeshape"} is a fully dichotomous binary tree. The purpose of the class \code{"treeshape"} is to study the topology of phylogenetic trees. The heights of branches are not provided for a tree of that class because we mainly focus on the balance aspect of the trees. The 'i'th row  of the nodes matrix represents the children of the node number i in the tree (\code{nodes[i,1]} being the left child, and \code{nodes[i,2]} being the right child). A positive value represents an internal node, while a negative one stands for a tip of the tree. The last row always represents the children of the root of the tree.
}

\value{
An object of class \code{"treeshape"}
}

\references{
Semple, C. and Steel, M. (2003) Phylogenetics. \emph{Oxford Lecture Series In Mathematics and its Applications}, \bold{24}, for the mathematical definitions and tools for phylogenetic trees.\cr
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}> \cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier François <\email{olivier.francois@imag.fr}> \cr
}

\seealso{
\code{\link{rtreeshape}}  
}

\examples{

## Nodes will define the nodes of a five tips tree
nodes<-matrix(nrow=4,ncol=2)
nodes[1,]<-c(-5,-4)
nodes[2,]<-c(1,-1)
nodes[3,]<-c(-3,2)
nodes[4,]<-c(-2,3)

## Now we can build the tree and plot it.
tree1<-treeshape(nodes)
plot(tree1)

## Computation of the sackin index for the tree :
sackin(tree1)

## Label will define the names of the tips
label=c("a", "b", "c", "d", "e")
tree2<-treeshape(nodes, label)
plot(tree1, tree2)
}
  
\keyword{datagen} 


