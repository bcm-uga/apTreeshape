\name{colless}
\alias{colless}

\title{ Compute the Colless' shape statistic on tree data }

\description{
\code{colless} computes the Colless' index of a tree and provides standardized values according to the Yule and PDA models.
}
\usage{
colless(tree, norm = NULL)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"} on which the Colless' index is computed.}
\item{norm}{ A character string equals to \code{NULL} (default) for no normalization or one of \code{"pda"} or \code{"yule"}.}
}

\details{
The Colless' index \eqn{Ic}{Ic} computes the sum of absolute values \eqn{|L-R|}{|L-R|} at each node of the tree where L (resp. R) is the size of the left (resp. right) daughter clade at the node.
  
The mean and standard deviation of the Colless's statistic \eqn{Ic}{Ic} have been computed by Blum et al (2005). Under the Yule model the standardized index \deqn{Iyule = \frac{Ic-n*\log(n)-n(\gamma-1-\log(2))}{n}}{Iyule = (Ic-n*log(n)-n(gamma-1-log(2)))/n} converges in distribution (\eqn{\gamma}{gamma} is the Euler constant). The limiting distribution is non Gaussian and is characterized as a functional fixed-point equation solution. Under the PDA model, the standardization is different \deqn{Ipda = \frac{Ic}{n^{3/2}}}{Ipda = Ic / n^(3/2)} and converges in distribution to the Airy distribution (See Flajolet and Louchard (2001)). Standardized indices are useful when one wishes to compare trees with different sizes. The \code{colless} function returns the value of the unnormalized index (default) or one of the standardized statistics (Yule or PDA).
  
}

\value{
An object of class \code{numeric} which is the Colless' index of the tree.
}

\references{ 
Mooers, A. O. and Heard, S. B. (1997) Inferring Evolutionnary Process from Phylogenetic Tree Shape. \emph{The Quarterly Review of Biology}, \bold{72}, 31 -- 54. 
  
Blum, M., Francois, O. and Janson, S. The mean, variance and limiting distribution of two statistics sensitive to phylogenetic tree balance; manuscript available from \cr \url{http://www-timc.imag.fr/Olivier.Francois/bfj.pdf}.

Flajolet, P. and Louchard, O. (2001) Analytic Variations on the Airy Distribution. \emph{Algorithmica}, \bold{31}, 361 -- 377.
}

\author{ 
Michael Blum <\email{michael.blum@imag.fr}>\cr
Nicolas Bortolussi <\email{nicolas.bortolussi@imag.fr}> \cr
Eric Durand <\email{eric.durand@imag.fr}>\cr
Olivier Francois <\email{olivier.francois@imag}>
}

\seealso{ \code{\link{sackin}} }

\examples{

## Colless' index for a randomly generated PDA tree (unnormalized value)
tpda<-rtreeshape(1,tip.number=70,model="pda")
colless(tpda[[1]],norm="pda")
  
## Histogram of Colless' indices for randomly generated Yule trees
main="Colless' indices for randomly generated Yule trees"
xlab="Colless' indices"
hist(sapply(rtreeshape(300,tip.number=50,model="yule"),FUN=colless,norm="yule"),
      freq=FALSE,main=main,xlab=xlab)
  
## Change the number of tips
hist(sapply(rtreeshape(300,tip.number=100,model="yule"),FUN=colless,norm="yule"),
      freq=FALSE,main=main,xlab=xlab)

}

\keyword{univar}
