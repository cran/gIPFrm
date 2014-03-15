\name{ipf.gamma}
\alias{ipf.gamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterative Proportional Fitting in Relational Models, with a Given Adjustment Factor
}
\description{
For a given model matrix and a given vector of observed cell frequencies, 
\code{ipf.gamma} computes the vector of frequencies whose subset sums are equal 
to the observed subset sums times the adjustment factor and whose  
relative frequencies satisfy the multiplicative structure prescribed by 
the model. 
}
\usage{
ipf.gamma(ModelMatrix, ObsTable, gamma, tol, estimand)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ModelMatrix }{an \code{I} by \code{J} model matrix of a relational model.
 \code{I} is the number of observations;  \code{J} is the number of generating subsets.  
}
  \item{ObsTable }{a vector of observed cell frequencies.
}
  \item{gamma }{an adjustment factor.
}
  \item{tol }{tolerance used in stopping criteria.
}
  \item{estimand }{set to \code{"probabilities"} in the case of multinomial 
      sampling; set to \code{"intensities"} in the case of Poisson sampling.
}
}


\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{model.matrix }{the model matrix.}
  \item{observed.data }{the vector of observed cell frequencies.}
  \item{fitted.values }{the estimated cell frequencies.}
  \item{model.parameters }{the estimated model parameters on the multiplicative scale.}
}
\references{
{ A.Klimova, T.Rudas, A.Dobra,
     Relational models for contingency tables. 
     J. Multivariate Anal., 2012, 104, 159--173.}

{A. Klimova, Coordinate-Free Exponential Families on Contingency Tables. PhD thesis. Advisers: Tamas Rudas and Thomas Richardson.}
}

\author{
 Anna Klimova, Tamas Rudas 
}




\examples{

## The model of independence for a 2 by 2 contingency table

A = matrix(c( 1,1,0,0,
              0,0,1,1,
              1,0,1,0,
              0,1,0,1), byrow=TRUE, nrow=4) ## the model matrix

y = c(1,2,3,4)  ## the observed data

ipf.gamma(A, y, 1, 1e-6, "intensities")   ## Adjustment factor is set to 1 

ipf.gamma(A, y, 0.5, 1e-6, "intensities")   ## Adjustment factor is set to 0.5 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Iterative Proportional Fitting }
