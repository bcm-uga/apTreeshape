\name{smaller.clade.spectrum}
\alias{smaller.clade.spectrum}

\title{ Compute the smaller clade spectrum of a tree. }

\description{
\code{smaller.clade.spectrum} returns a n*2 matrix where n is the number of internal nodes of the tree. For each i in 1:n, the [i,1] element of the matrix is the size of the clade rooted at the 'i'th node of the tree. [i,2] is the size of the smaller daughter clade of the 'i'th node of the tree.
}

\usage{
smaller.clade.spectrum(tree)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
}

\value{
A n*2 matrix (where n is the number of internal nodes of the tree) containing the size of the clades and the smaller clades.
\code{smaller.clade.spectrum(tree)[1,1]} contains the number of tips of the tree.
\code{smaller.clade.spectrum(tree)[i,1]} contains the number of tips of the subtree whose root is the node number \code{n-i+1}.
\code{smaller.clade.spectrum(tree)[1,2]} contains the number of tips of the smaller daughter clade at the root.  
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.franois@imag.fr}>
}

\seealso{
\code{\link{spectrum.treeshape}} 
}

\examples{

# computes the log-likelihood for Aldous' model
shape.new <- function(tree){
h <- function(n){sum(1/(1:n))}
mat <- smaller.clade.spectrum(tree)
parent <- mat[,1]
daughter <- mat[,2]
nh.n <- sapply(parent, FUN = function(x){x*h(x-1)} )
s <- sum(log(daughter/parent) + log(1 - daughter/parent) + log(nh.n))
return(s)}
 
# distribution over 200 replicates

tr <- rtreeshape(200, 100, FUN =
function(n,i){if((i>0)&(i<n))return(1/i/(n-i)) else return(0)})
res <- sapply( tr, FUN = shape.new)  
hist(res)}

\keyword{univar}