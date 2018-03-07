\name{lambda_N}
\alias{lambda_N}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computation of lambda_epsilon
}
\description{
Computes the parameter lambda_epsilon from the online appendix for n in 1:N
}
\usage{
lambda_N(epsilon, beta, N, upper = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epsilon}{
Minimum size of unsampled splits (see appendix 1)
}
  \item{beta}{
The beta parameter quantifying the balance of the tree
}
  \item{N}{
Tip number
}
  \item{upper}{
Upper bound from the integral (see lambda_n expression, default to 10000)
}
}
\author{
Odile Maliet, Fanny Gascuel & Amaury Lambert
}


