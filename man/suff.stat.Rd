\name{suff.stat}
\alias{suff.stat}
\title{
 Sufficient Statistics under a Relational Model 
}
\description{
  For an \code{I} by \code{J} model matrix of a relational model and 
  a vector of frequencies of the length \code{I}, the function computes 
  sufficient statistics under the model (subset sums). 
}


\usage{
suff.stat(ModelMatrix, Table)
}
\arguments{
  \item{ModelMatrix}{
     \code{ModelMatrix} a model matrix of a relational model. 
}
  \item{Table}{
      \code{Table} a vector of frequencies. 
}
}

\value{ a vector of subset sums.}

\references{
     A.Klimova, T.Rudas, A.Dobra,
     Relational models for contingency tables. 
     J. Multivariate Anal., 104, 159--173. 
}
\author{
  Anna Klimova
}



\examples{

### Multiplicative model from Aitchison and Silvey (1960)

A = matrix(c(1, 0, 0, 1, 0, 1, 1, 
             0, 1, 0, 1, 1, 0, 1,
             0, 0, 1, 0, 1, 1, 1), byrow=TRUE, nrow=3) ## the model matrix 

y = c(46,24,7,15,3,4,1) ## the observed data

suff.stat(A,y) 

### pneumonia infection in calves:

A=matrix(c(2,1,0,0,1,1), 2,3,byrow=TRUE) ## the model matrix

y=c(30,63,63)  ## the observed data

suff.stat(A,y) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Sufficient Statistics }
