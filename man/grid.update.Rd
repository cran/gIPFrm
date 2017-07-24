\name{grid.update}
\alias{grid.update}
\title{
 Search for the adjustment factor corresponding to the MLE using a grid
}
\description{
\code{grid.update} computes the segment 
that is known to contain the adjustment factor corresponding to the MLE
and, using a grid on this segment, finds this adjsutment factor. It is needed only for relational models for probabilities.
}
\usage{
grid.update(ModelMx, ObsTbl, tolerance)
}


\arguments{
  \item{ModelMx }{an \code{I} by \code{J} model matrix of a relational model. Here
 \code{I} is the number of observations and \code{J} is the number of generating subsets.  
}
  \item{ObsTbl }{a vector of observed cell frequencies of length \code{I}. 
}
  \item{tolerance }{tolerance used in stopping criteria.
}
}

\value{
  \item{gamma.tilde }{the adjustment factor under the precision given by \code{tolerance}.}
  \item{model.tilde }{the value returned by \code{ipf.gamma()}
                      with the adjustment factor \code{gamma} equal to \code{gamma.tilde}.}
}

\references{

{A. Klimova, Coordinate-Free Exponential Families on Contingency Tables. PhD thesis. Advisers: Tamas Rudas and Thomas Richardson.}

}

\author{ Anna Klimova, Tamas Rudas}

\examples{

### Multiplicative model from Aitchison and Silvey (1960)

A = matrix(c(1, 0, 0, 1, 0, 1, 1, 
             0, 1, 0, 1, 1, 0, 1,
             0, 0, 1, 0, 1, 1, 1), byrow=TRUE, nrow=3) ## the model matrix 

y = c(46,24,7,15,3,4,1) ## the observed data

grid.update(A, y, 1e-4)


## The model of independence for a 2 by 2 contingency table

A = matrix(c( 1,1,0,0,
              0,0,1,1,
              1,0,1,0,
              0,1,0,1), byrow=TRUE, nrow=4) ## the model matrix

y = c(1,2,3,4)  ## the observed data


grid.update(A, y, 1e-5)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
