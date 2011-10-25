\name{skmeans_xdist}
\alias{skmeans_xdist}
\title{Cosine Cross-Distances}
\description{
  Compute cosine cross-distances between the rows of matrices.
}
\usage{
skmeans_xdist(x, y = NULL)
}
\arguments{
  \item{x}{A numeric data matrix.
    Can be a dense matrix,
    \link[slam:matrix]{simple triplet matrix}
    (package \pkg{slam}), or a
    \link[Matrix:dgTMatrix-class]{dgTMatrix}
    (package \pkg{Matrix}).}
  \item{y}{\code{NULL} (default), or as for \code{x}.
    The default is equivalent to taking \code{y} as \code{x} (but more
    efficient).}
}
\value{
  A dense matrix \eqn{d} with entry \eqn{d_{ij} = 1 - \cos(x_i, y_j)}
  the cosine distance between the \eqn{i}-th row \eqn{x_i} of \code{x}
  and the \eqn{j}-th row \eqn{y_j} of \code{y}.
}
