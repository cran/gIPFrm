\name{single.cells}
\alias{single.cells}
\title{
Single Cells under a Relational Model
}
\description{
The function finds all single cells under a relational model. Such cells appear 
as the only positive entries in their row and column in the model matrix. 
}
\usage{
single.cells(ModelMatrix)
}

\arguments{
  \item{ModelMatrix }{a model matrix of a relational model.}
}

\value{
{the row and column indices of the single cells.}
}

\author{
  Anna Klimova
}


\examples{
G = matrix(c( 2,0,0,0,0,
              0,1,1,0,0,
              0,0,1,1,1), byrow=TRUE, nrow=3) ## a given matrix

single.cells(G)  

}

