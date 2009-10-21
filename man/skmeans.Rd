\name{skmeans}
\alias{skmeans}
\title{Compute Spherical k-Means Partitions}
\description{
  Partition given vectors \eqn{x_b} by minimizing the criterion
  \eqn{\sum_{b,j} w_b u_{bj}^m d(x_b, p_j)}
  where the \eqn{w_b} are case weights,
  \eqn{u_{bj}} is the membership of \eqn{x_b} to class \eqn{j},
  \eqn{p_j} is the \emph{prototype} of class \eqn{j}
  (thus minimizing \eqn{\sum_b w_b u_{bj}^m d(x_b, p)} over \eqn{p}),
  and \eqn{d} is the cosine dissimilarity
  \eqn{d(x, p) = 1 - \cos(x, p)}.
}
\usage{
skmeans(x, k, method = NULL, m = 1, weights = 1, control = list())
}
\arguments{
  \item{x}{A numeric data matrix, with rows corresponding to the objects
    to be partitioned (such that row \eqn{b} contains \eqn{x_b}).  Can
    be a dense matrix, a
    \link[slam:matrix]{simple triplet matrix}
    (package \pkg{slam}), or a
    \link[Matrix:dgTMatrix-class]{dgTMatrix}
    (package \pkg{Matrix}).}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
  \item{method}{a character string specifying one of the built-in
    methods for computing spherical \eqn{k}-means partitions, or a
    function to be taken as a user-defined method, or \code{NULL}
    (default value).  If a character string, its lower-cased version is
    matched against the lower-cased names of the available built-in
    methods using \code{\link{pmatch}}.  See \bold{Details} for
    available built-in methods and defaults.}
  \item{m}{a number not less than 1 controlling the softness of the
    partition (as the \dQuote{fuzzification parameter} of the fuzzy
    \eqn{c}-means algorithm).  The default value of 1 corresponds to
    hard partitions obtained from the standard spherical \eqn{k}-means
    problem; values greater than one give partitions of increasing
    softness obtained from a generalized fuzzy spherical \eqn{k}-means
    problem.}
  \item{weights}{a numeric vector of non-negative case weights.
    Recycled to the number of elements given by \code{x} if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\details{
  The \dQuote{simple} spherical \eqn{k}-means problem where all case
  weights are one and \eqn{m = 1} is equivalent to maximizing the
  criterion \eqn{\sum_j \sum_{b \in C_j} \cos(x_b, p_j)},
  where \eqn{C_j} is the \eqn{j}-th class of the partition.  This is the
  formulation used in Dhillon & Modha (2001) and related references, and
  is criterion function \eqn{\mathcal{I}_2} in the CLUTO documentation.

  Obtaining optimal spherical \eqn{k}-means partitions obviously is a
  computationally hard problem, and several methods are available which
  attempt to obtain optimal partitions.  The built-in methods are as
  follows.

  \describe{
    \item{\code{"genetic"}}{a genetic algorithm patterned after the
      genetic \eqn{k}-means algorithm of Krishna & Narasimha Murty (1999).}
    \item{\code{"pclust"}}{a Lloyd-Forgy style fixed-point algorithm which
      iterates between determining optimal memberships for fixed
      prototypes, and computing optimal prototypes for fixed
      memberships, using the general-purpose prototype-based
      partitioning framework of \code{\link[clue]{pclust}} in package
      \pkg{clue}.}
    \item{\code{"CLUTO"}}{an interface to the \code{vcluster} partitional
      clustering program from CLUTO, the CLUstering TOolkit by George
      Karypis.}
    \item{\code{"lih"}}{the Local Improvement Heuristic of Dhillon, Guan
      and Kogan (2002).  If a fixed-point iteration does not (sufficiently)
      improve the criterion, improvement by swapping the clusters of two
      elements is attempted.}
    \item{\code{"lihc"}}{the Local Improvement Heuristic with Chains
      method of Dhillon, Guan and Kogan (2002).  This attempts further
      improvements via Kernighan-Lin object moves chains.}
  }

  Method \code{"pclust"} is the only method available for soft spherical
  \eqn{k}-means problems.  Method \code{"genetic"} can handle case
  weights.  By default, the genetic algorithm is used for obtained hard
  partitions, and the fixed-point algorithm otherwise.

  Control parameters for method \code{"genetic"} are as follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximum number of
      iterations for the genetic algorithm.  Defaults to 12.}
    \item{\code{popsize}}{an integer giving the population size for the
      genetic algorithm.  Default: 6.}
    \item{\code{mutations}}{a number between 0 and 1 giving the
      probability of mutation per iteration.  Defaults to 0.1.}
    \item{\code{start}}{a list with the prototypes for the initial
      population.}
    \item{\code{reltol}}{The minimum relative improvement per
      iteration.  If improvement is less, the algorithm will stop under
      the assumption that no further significant improvement can be
      made.  Defaults to 1e-8.}
    \item{\code{verbose}}{a logical indicating whether to provide
      some output on minimization progress.
      Defaults to \code{getOption("verbose")}.}
  }

  See the documentation of \code{\link[clue]{pclust}} in package
  \pkg{clue} for the control parameters for method \code{"pclust"}.

  Control parameters for method \code{"CLUTO"} are as follows.
  \describe{
    \item{\code{vcluster}}{the path to the CLUTO \code{vcluster}
      executable.}
    \item{\code{colmodel}}{a specification of the CLUTO column model.
      See the CLUTO documentation for more details.}
    \item{\code{verbose}}{as for the genetic algorithm.}
    \item{\code{control}}{a character string specifying arguments passed
      over to the \code{vcluster} executable.}
  }

  Control parameters for the local improvement heuristics are as
  follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be performed.
      Defaults to 100.}
    \item{\code{start}}{a single prototype to be used as a starting 
      value.}
    \item{\code{reltol}}{as for the genetic algorithm.}
    \item{\code{verbose}}{as for the genetic algorithm.}
  }
  The enhanced chain-based heuristic has the following additional
  control parameter.
  \describe{
    \item{\code{maxchains}}{an integer giving the maximal length of the
      Kernighan-Lin chains.  Defaults to 10.}
  }
  
  Method \code{"CLUTO"} requires that the CLUTO \code{vcluster}
  executable is available.  CLUTO binaries for the Linux, Sun, OSX, and
  MS Windows platforms can be obtained from
  \url{http://www-users.cs.umn.edu/~karypis/cluto/}.

  User-defined methods must have formals \code{x}, \code{k} and
  \code{control}, and optionally may have formals \code{weights} 
  or \code{m} if providing support for case weights or soft spherical
  \eqn{k}-means partitions, respectively.
}
\value{
  An object of class \code{pclust} (see the information on
  \link[clue:pclust]{pclust objects} in package \pkg{clue} for
  further details) representing the obtained spherical \eqn{k}-means
  partition, which is a list with components including the following:
  \item{prototypes}{a dense matrix with \code{k} rows giving the
    prototypes.}
  \item{membership}{cluster membership as a matrix with \code{k}
    columns.}
  \item{cluster}{the class ids of the closest hard partition (the
    partition itself if \eqn{m = 1}).}
  \item{value}{the value of the criterion.}
}
\references{
  I. S. Dhillon and D. S. Modha (2001).
  Concept decompositions for large sparse text data using clustering.
  \emph{Machine Learning}, \bold{42}, 143--175.

  I. S. Dhillon and Y. Guan and J. Kogan (2002).
  Iterative clustering of high dimensional text data augmented by local
  search.
  In \emph{Proceedings of the Second IEEE International Conference on
    Data Mining}, pages 131--138.
  \url{http://www.cs.utexas.edu/users/inderjit/public_papers/iterative_icdm02.pdf}.

  K. Krishna and M. Narasimha Murty (1999).
  Genetic \eqn{K}-means algorithm.
  \emph{IEEE Transactions on Systems, Man, and Cybernetics --- Part B:
    Cybernetics}, \bold{29}/3, 433--439.
  \url{http://eprints.iisc.ernet.in/2937/1/genetic-k.pdf}.

  G. Karypis (2003).
  \emph{CLUTO: A Clustering Toolkit}.
  Technical Report #02-017, Department of Computer Science, University
  of Minnesota.
  \url{http://glaros.dtc.umn.edu/gkhome/fetch/sw/cluto/manual.pdf}.
}
\author{
  Kurt Hornik \email{Kurt.Hornik@wu.ac.at}, \cr
  Ingo Feinerer \email{feinerer@logic.at}, \cr
  Martin Kober \email{martin.kober@wu.ac.at}.
}
\examples{
## Use CLUTO dataset 're0':
x <- readCM(system.file("cluto", "re0.mat",
                        package = "skmeans"),
            system.file("cluto", "re0.mat.clabel",
                        package = "skmeans"))
## Which is not really small:
dim(x)
## Partition into 5 clusters.
party <- skmeans(x, 5, control = list(verbose = TRUE))
## Criterion value obtained:
party$value
## Compare with "true" classifications:
class_ids <-
    readLines(system.file("cluto", "re0.mat.rclass",
                          package = "skmeans"))
table(class_ids, party$cluster)
}
\keyword{cluster}
