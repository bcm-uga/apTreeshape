\name{get_PD_sample}
\alias{get_PD_sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the proportion of conserved PD
}
\description{
Computes the proportion of conserved phylogenetic diversity as a function of the proportion of conserved species, in trees simulated by the model
}
\usage{
get_PD_sample(epsilon, beta, alpha, N, sampl.frac, ntree, equal.ab,
              eta, lengths = "yule", b = 1, d = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
Imbalance index
}
  \item{alpha}{
Clade age-richness index
}
  \item{N}{
Initial tip number
}
  \item{sampl.frac}{
Vector of tips fractions for which we want to compute the conserved PD
}
  \item{ntree}{
Number of simulated trees
}
  \item{equal.ab}{
If set to TRUE, all species have the same probability to go extinct first (default to TRUE)
}
  \item{eta}{
Clade abundance-richness index (if equal.ab == FALSE)
}
  \item{lengths}{
Model used to simulate node depths (can be "yule" (the default) or "kingman")}
  \item{b}{
Birth rate (if lengths == "yule")
}
  \item{d}{
Death rate (if lengths == "yule")
}
}
\details{
A table of size length(sample.fract)*ntree. The element in ligne i and column j is the remaining Phylogenetic Diversity fraction for the j^th tree in wich a fraction sample.frac[i] has been sampled.  
}
\references{
Maliet O., Gascuel F., Lambert A. (2018) \emph{Ranked tree shapes, non-random 
extinctions and the loss of phylogenetic diversity}, bioRxiv 224295, 
doi: https://doi.org/10.1101/224295 
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}
\examples{
set.seed(813)
sampl.frac=seq(0,1,by=0.01) 

PD = get_PD_sample(epsilon=0.01,beta=0,alpha=0,N=100,
                sampl.frac=sampl.frac,ntree=100,equal.ab=FALSE,
                eta=0.5,lengths="yule",b=1,d=0)
probs = c(0.025, 0.975, 0.5)  
PD_stats = as.vector(t(sapply(1:nrow(PD), function(i){quantile(PD[i,], probs)}))) 
par(mgp=c(2.2, 0.8, 0))
par(mar=c(4, 4, 1, 1))
plot(1, type="n", xlab="Fraction of extinct species, p",
    ylab="PD loss", ylim=c(0,1), xlim=c(0,1))
points(c(0,1),c(0,1),t="l",col=grey(0),lty=3)
# plot 95% confidence intervals (grey area)
polygon(c(1-sampl.frac, rev(1-sampl.frac)), c(1-PD_stats[(1:length(sampl.frac))],
      rev(1-PD_stats[((length(sampl.frac)+1):(2*length(sampl.frac)))])),
      border=NA, col=grey(0.7))
# plot median value (black line) 
points(1-sampl.frac, 1-PD_stats[((2*length(sampl.frac)+1):(3*length(sampl.frac)))],t="l")

}
