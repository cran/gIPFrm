\name{gIPFrm-package}
\alias{gIPFrm-package}
\alias{gIPFrm}
\docType{package}
\title{
Generalized Iterative Proportional Fitting for Relational Models 
}
\description{
 The package provides an iterative scaling procedure that computes the maximum likelihood estimates of the cell frequencies and of the model
 parameters under a relational model, with or without the overall effect. }
\details{
\tabular{ll}{
Package: \tab gIPFrm\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-03-21\cr
License: \tab GPL (>= 2)\cr
}
The iterative proportional fitting procedure is called by the function 
\code{g.ipf}. 
}
\author{
Anna Klimova, Tamas Rudas

Maintainer: 
Anna Klimova <aklimova25@gmail.com> 
}
\references{
{A.Klimova, T.Rudas, A.Dobra,
     Relational models for contingency tables. 
     J. Multivariate Anal., 104, 159--173.} 

{A. Klimova, Coordinate-Free Exponential Families on Contingency Tables. PhD thesis. Advisers: Tamas Rudas and Thomas Richardson.}

{A.Agresti, Categorical Data Analysis. Wiley, New York, 1990.}

{J.Aitchison, S.D.Silvey, Maximum-likelihood estimation procedures and associated tests of significance.
 J. Roy. Statist. Soc. Ser.B, 1960, 22, 154--171.}

{G.Kawamura, T.Matsuoka, T.Tajiri, M.Nishida, M.Hayashi,
 Effectiveness of a sugarcane-fish combination as bait in trapping swimming crabs.
 Fisheries Research, 1995, 22, 155--160.}

}

\note{
Tamas Rudas was supported in part by Grant K-106154 
from the Hungarian National Scientific Research Fund (OTKA).
The authors wish to thank Juraj Medzihorsky for his help 
with building this package.
}

\examples{

### Multiplicative model from Aitchison and Silvey (1960)

A = matrix(c(1, 0, 0, 1, 0, 1, 1, 
             0, 1, 0, 1, 1, 0, 1,
             0, 0, 1, 0, 1, 1, 1), byrow=TRUE, nrow=3) ## the model matrix 

y = c(46,24,7,15,3,4,1) ## the observed data

g.ipf(A, y, 1e-4, "probabilities", "grid")
g.ipf(A, y, 1e-6, "probabilities", "bisection")

}
