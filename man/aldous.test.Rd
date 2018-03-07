\name{aldous.test}
\alias{aldous.test}

\title{ Visualizing balance via scatter diagrams }

\description{
A graphical test to decide if tree data fit the Yule or the PDA  models.
}

\usage{
aldous.test(tree, xmin=20, \dots)
}

\arguments{
\item{tree}{ An object of class \code{"treeshape"}.}
\item{xmin}{ An object of class \code{"numeric"} that defines the range of the x-axis. The minimal parent clade size displayed in the graphical representation (default: \code{xmin=20}).}
\item{\dots}{ further arguments to be passed to \code{plot()}.}
}

\details{
A binary tree contains a set of splits \deqn{(m,i) = (size~of~parent~clade, size~of~smaller~daughter~clade)}{(m,i)=(size of parent clade, size of smaller daughter clade)} which can be plotted as a scatter diagram. Aldous' proposal for studying tree balance is that, given a large phylogenetic tree, one should estimate the median size of the smaller daughter clade as a function of the parent clade and use this function as a descriptor of balance or imbalance of the tree. It is convenient to make a log-log plot and to ignore small parent clades. The scatter diagram shows lines giving the approximate median values of the size of smaller daughter clade predicted by the beta-splitting model for two values of beta, the value for the Yule \eqn{(\beta=0)}{(beta=0)} and PDA (\eqn{\beta=-1.5}{beta=-1.5}) models. In other words, if the null model were true, then the scatter diagram for a typical tree would have about half the points above the line and half below the line, throughout the range.\cr
  
The green line represents the median regression estimated from the tree data.
}

\value{
The function provides a graphical display of results.
}

\references{ 
Aldous, D. J. (1996) \emph{Probability Distributions on Cladograms}. pp.1-18 of \emph{Random Discrete Structures} eds D. Aldous and R. Pemantle, IMA Volumes Math. Appl. 76. 
  
Aldous, D. J. (2001) Stochastic Models and Descriptive Statistics for Phylogenetic Trees, from Yule to Today. \emph{Statistical Science}, \bold{16}, 23 -- 24. 
}

\author{ 
Michael Blum \email{michael.blum@imag.fr}\cr
Nicolas Bortolussi \email{nicolas.bortolussi@imag.fr}\cr
Eric Durand \email{eric.durand@imag.fr}\cr
Olivier Francois \email{olivier.francois@imag.fr}
}

\examples{
library(quantreg)
aldous.test(rbiased(2000, p=.5))

## Test with a huge balanced tree:
aldous.test(rbiased(2000, p=.5))


  
}

\keyword{htest}
