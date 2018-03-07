\name{sackin}
\alias{sackin}

\title{ Compute the Sackin's index of a tree }

\description{
\code{sackin} computes the Sackin's index on tree and normalizes it.
}
\usage{
sackin(tree, norm = NULL)
}

\arguments{
\item{tree}{An object of class \code{"treeshape"} on which the Sackin's index is computed. }
\item{norm}{A character string equals to \code{"null"} (default) for no normalization, \code{"pda"} for the PDA normalization or \code{"yule"} for the Yule normalization. }
}

\details{
The Sackin's index is computed as the sum of the number of ancestors for each tips of the tree. The less balanced a tree is and the larger its Sackin's index. It can be normalized in order to obtain a statistic that does not depend on the tree size, and so compare  trees with different sizes. The normalization depends on the reference model (Yule or PDA). Under the Yule model, the normalized index is \deqn{Iyule = \frac{Is-2n*\sum_{j=2}^n{\frac{1}{j}}}{n}}{Iyule = (Is-2*n*sum1/j,j=2..n))/n} where \eqn{Is}{Is} is the non-normalized Sackin's index for a n-tips tree. Under the PDA model, the normalized index is \deqn{Ipda = \frac{Is}{n^{3/2}}}{Ipda = Is/ n^(3/2)} See details on the Colless index.
}

\value{
An object of class \code{numeric} which contains the Sackin's index of the tree.
}

\references{
Mooers, A. O., Heard, S. B. (Mar., 1997) Inferring Evolutionnary Process from Phylogenetic Tree Shape. \emph{The Quarterly Review of Biology}, \bold{72}, 31 -- 54, for more details about the Sackin'index and its signification about the balance of phylogenetic trees. 
  
Blum, M., Francois, O. and Janson, S. The mean, variance and limiting distribution of two statistics sensitive to phylogenetic tree balance; manuscript available from \cr \url{http://www-timc.imag.fr/Olivier.Francois/bfj.pdf}.
}

\author{
Michael Blum <\email{michael.blum@imag.fr}> \cr 
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}> \cr
Olivier François <\email{olivier.francois@imag.fr}>
}

\seealso{
\code{\link{colless}}\cr
\code{\link{sackin.test}}  
}

\examples{

## Index of Sackin of a PDA tree :
tpda<-rtreeshape(1,tip.number=70,model="pda")
tpda<-tpda[[1]]
sackin(tpda,norm="pda")
  
## Histogram of the Sackin's indices for randomly generated Yule trees, 
##      with no normalization
main="Histogram of Sackin's indices for randomly generated Yule trees"
xlab="Sackin's index"
hist(sapply(rtreeshape(300,tip.number=50,model="yule"),FUN=sackin,norm="yule"),
      freq=FALSE, main=main, xlab=xlab)

## Change the size of the trees:
hist(sapply(rtreeshape(300,tip.number=100,model="yule"),FUN=sackin,norm="yule"),
      freq=FALSE, main=main, xlab=xlab)
}

\keyword{univar}
