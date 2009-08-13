\name{cluto}
\alias{readCM}
\alias{writeCM}
\title{Read and Write Matrices in CLUTO Sparse Matrix Format}
\description{
  Read and write matrices in CLUTO sparse matrix format.
}
\usage{
readCM(con, clabel = NULL)
writeCM(mat, con)
}
\arguments{
  \item{con}{a connection object or a character string.}
  \item{clabel}{a character string giving the path to a CLUTO column
    labels file, or \code{NULL} (default).}
  \item{mat}{a matrix object.}
}
\details{
  Documentation for CLUTO including its sparse matrix format is
  available from \url{http://www-users.cs.umn.edu/~karypis/cluto/}.

  \code{readCM} reads CLUTO sparse matrices, returning a simple triplet
  matrix (see \code{\link[slam]{simple_triplet_matrix}} in package
  \pkg{slam} for details).

  \code{writeCM} writes CLUTO sparse matrices.  Argument \code{mat} must
  be coercible to a simple triplet matrix via
  \code{\link[slam]{as.simple_triplet_matrix}}.
}
\keyword{IO}
