\name{g.ipf}
\alias{g.ipf}
\title{
 Generalized Iterative Proportional Fitting for Relational Models
}
\description{
\code{g.ipf} computes the maximum likelihood estimates of the cell frequencies and of the model 
parameters under a relational model specified by a model matrix.  
}
\usage{
g.ipf(ModelMatrix, ObsTable, tol, estimand, adjustment)
}


\arguments{
  \item{ModelMatrix }{an \code{I} by \code{J} model matrix of a relational model. Here
 \code{I} is the number of observations and \code{J} is the number of generating subsets.  
}
  \item{ObsTable }{a vector of observed cell frequencies of length \code{I}. 
}
  \item{tol }{tolerance used in stopping criteria.
}
  \item{estimand }{set to \code{"probabilities"} in the case of multinomial 
      sampling; set to \code{"intensities"} in the case of Poisson sampling.   
}
  \item{adjustment }{set to \code{"grid"} if a grid is 
      used to update the adjustment factor; set to \code{"bisection"} if the bisection method is used to update the 
      adjustment factor; 
      set to \code{"none"} if \code{estimand} is \code{"intensities"}.
}
}

\value{
  \item{model.matrix }{the model matrix.}
  \item{observed.data }{the vector of observed cell frequencies.}
  \item{fitted.values }{the maximum likelihood estimates of the cell frequencies.}
  \item{estimated.total }{the sum of the estimated cell frequencies.}
  \item{adjustment.for.total }{the estimated total divided by the observed total.}
  \item{adjustment.for.subsets }{the adjustment factor for the subset sums.}
  \item{model.parameters }{the maximum likelihood estimates for model parameters on the multiplicative scale.}
  \item{degrees.of.freedom }{the degrees of freedom of the relational model.}
  \item{chisq.statistic }{Pearson's chi-squared statistic.}
  \item{p.value.chisq }{the p-value, based on Pearson's chi-squared statistic.}
  \item{log.likelihood.ratio.statistic }{the log likelihood ratio statistic.}
  \item{p.value.log.likelihood.ratio.statistic }{the p-value, based on the log likelihood ratio statistic.}
  \item{Bregman.statistic }{the Bregman statistic.}
  \item{p.value.Bregman.statistic }{the p-value, based on the Bregman statistic.}
}

\references{
{A.Klimova, T.Rudas, A.Dobra,
     Relational models for contingency tables. 
     J. Multivariate Anal., 2012, 104, 159--173.} 

{A.Klimova, T.Rudas, 
     Iterative proportional scaling for curved exponential families. 
     Scand. J. Statist., 2015, 42, 832--847.} 

{A. Klimova, Coordinate-Free Exponential Families on Contingency Tables. PhD thesis. Advisers: Tamas Rudas and Thomas Richardson.}

{A.Agresti, Categorical Data Analysis. Wiley, New York, 1990.}

{J.Aitchison, S.D.Silvey, Maximum-likelihood estimation procedures and associated tests of significance.
 J. Roy. Statist. Soc. Ser.B, 1960, 22, 154--171.}

{G.Kawamura, T.Matsuoka, T.Tajiri, M.Nishida, M.Hayashi,
 Effectiveness of a sugarcane-fish combination as bait in trapping swimming crabs.
 Fisheries Research, 1995, 22, 155--160.}

}

\author{ Anna Klimova, Tamas Rudas}


\examples{

### Multiplicative model from Aitchison and Silvey (1960)

A = matrix(c(1, 0, 0, 1, 0, 1, 1, 
             0, 1, 0, 1, 1, 0, 1,
             0, 0, 1, 0, 1, 1, 1), byrow=TRUE, nrow=3) ## the model matrix 

y = c(46,24,7,15,3,4,1) ## the observed data

g.ipf(A, y, 1e-6, "probabilities", "bisection")
g.ipf(A, y, 1e-4, "probabilities", "grid")



### Bait study in swimming crabs, see Kawamura at al. (1995) :  

A = matrix(c(1,1,0,1,0,1), 2,3,byrow=TRUE) ## the model matrix

y1 = c(36,2,11)  ## the observed data for Charybdis japonica
y2 = c(71,3,44)  ## the observed data for Portunuspelagicus 

## If the model is for intensities, the adjustemnt factor is equal to one;
## the adjustment parameter is redundant.

g.ipf(A, y1, 1e-6, "intensities", "none")
g.ipf(A, y2, 1e-6, "intensities", "none")



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Exponential Families }
\keyword{ Iterative Proportional Fitting }% __ONLY ONE__ keyword per line
