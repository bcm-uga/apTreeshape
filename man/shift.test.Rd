\name{shift.test}
\alias{shift.test}

\title{ Testing diversification rate variation in phylogenetic trees }

\description{
A statistical test of diversification rate shift within an evolutionary process.
}

\usage{
shift.test(tree, node, lambda1, lambda2, nrep, silent)
}

\arguments{
\item{tree}{An object of class \code{treeshape}.}
\item{node}{An integer value that represents the internal node at which the test is done. It ranges from 1 to the number of tips minus 2.}
\item{lambda1}{A positive numerical value (object of class \code{numeric}) that represents an ancestral diversification rate in the lineage ending at \code{node}.}
\item{lambda2}{A positive numerical value (object of class \code{numeric}) that represents the shifted diversification rate after the speciation event at \code{node}.}
\item{nrep}{Number of Monte-carlo replicates in the computation of the P-value.} 
\item{silent}{A boolean indicating whether the test is silent or verbose.}
}

\details{This function implements a test of diversification rate shift based on the Delta1 statistic (Moore et al. 2004). We introduced some simplifications of the test statistic
using basic results in branching process.}

\value{
The function provides textual results and a P-value for the null hypothesis of no diversification rate shift against a (lambda2/lambda1)-fold shift at the node under consideration.
}

\references{ 
Moore, B.R., Chan, K.M.A, and Donoghue M.J. (2004) \emph{Detecting diversification rate variation in supertrees}. In Bininda-Emonds, O. R. P. (ed.) \emph{Phylogenetic Supertrees: Combining Information to Reveal the Tree of Life}, pp. 487-533. Computational Biology, volume 3 (Dress, A., series ed.).
}

\author{ 
Eric Durand \email{eric.durand@imag.fr}\cr
Olivier Francois \email{olivier.francois@imag.fr}
}

\examples{
	## Detecting diversification rate variation in bird families (135 tips)
 data(bird.families)
 tree.birds <- as.treeshape(bird.families, model = "yule")
 class(tree.birds) <- "treeshape"
 pv <- sapply(1:135, FUN = function(i) shift.test(tree.birds, i, lambda1 = 1, 
                                      lambda2 = 100, nrep = 1000, silent = TRUE))

	## Significant shifts detected at nodes = 67 and 78	
 pv[c(67,78)]
 shift.test(tree.birds, node = 67, lambda1 = 1, lambda2 = 100, nrep = 10000, silent = TRUE)
 shift.test(tree.birds, node = 78, lambda1 = 1, lambda2 = 100, nrep = 10000, silent = TRUE)

 	## visualize the shifts
 par(mfrow=c(2,1)) 
 plot(cutreeshape(tree.birds, ancestor(tree.birds, 67) , "bottom"))
 plot(cutreeshape(tree.birds, 78 , "bottom"))

}

\keyword{htest}
