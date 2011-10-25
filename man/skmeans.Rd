\name{skmeans}
\alias{skmeans}
\title{Compute Spherical k-Means Partitions}
\description{
  Partition given vectors \eqn{x_b} by minimizing the spherical
  \eqn{k}-means criterion
  \eqn{\sum_{b,j} w_b u_{bj}^m d(x_b, p_j)}
  over memberships and prototypes,
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
    (package \pkg{Matrix}).
    Zero rows are not allowed.
  }
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
    hard partitions; values greater than one give partitions of
    increasing softness obtained from a generalized soft spherical
    \eqn{k}-means problem.}
  \item{weights}{a numeric vector of non-negative case weights.
    Recycled to the number of objects given by \code{x} if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\details{
  The \dQuote{standard} spherical \eqn{k}-means problem where all case
  weights are one and \eqn{m = 1} is equivalent to maximizing the
  criterion \eqn{\sum_j \sum_{b \in C_j} \cos(x_b, p_j)},
  where \eqn{C_j} is the \eqn{j}-th class of the partition.  This is the
  formulation used in Dhillon & Modha (2001) and related references, and
  when optimized over the prototypes yields the criterion function
  \eqn{\mathcal{I}_2} in the CLUTO documentation.

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
      memberships.  For hard partitions, this can optionally attempt
      further local improvements via Kernighan-Lin chains of first
      variation single object moves as suggested by Dhillon, Guan and
      Kogan (2002).}
    \item{\code{"CLUTO"}}{an interface to the \code{vcluster} partitional
      clustering program from CLUTO, the CLUstering TOolkit by George
      Karypis.}
    \item{\code{"gmeans"}}{an interface to the \code{gmeans} partitional
      clustering program by Yuqiang Guan.}
    \item{\code{"kmndirs"}}{an interface to the C code for the
      \eqn{k}-mean-directions algorithm of Ranjan Maitra and Ivan
      P. Ramler.}
  }

  Method \code{"pclust"} is the only method available for soft spherical
  \eqn{k}-means problems.  Method \code{"genetic"} can handle case
  weights.  By default, the genetic algorithm is used for obtaining hard
  partitions, and the fixed-point algorithm otherwise.

  Common control parameters for methods \code{"genetic"} and
  \code{"pclust"} are as follows.
  \describe{
    \item{\code{start}}{a specification of the starting values to be
      employed.  Can either be a character vector with elements
      \code{"p"} (randomly pick objects as prototypes), \code{"i"}
      (randomly pick ids for the objects), \code{"S"} (take \eqn{p}
      minimizing \eqn{\sum_b w_b d(x_b, p)} as the first prototype, and
      successively pick objects farthest away from the already picked
      prototypes), or \code{"s"} (like \code{"S"}, but with the first
      prototype a randomly picked object).  Can also be a list of
      \code{skmeans} objects (obtained by previous runs), a list of
      prototype matrices, or a list of class ids.  For the genetic
      algorithm, the given starting values are used as the initial
      population; the fixed-point algorithm is applied individually to
      each starting value, and the best solution found is returned.}
    \item{\code{reltol}}{The minimum relative improvement per
      iteration.  If improvement is less, the algorithm will stop under
      the assumption that no further significant improvement can be
      made.  Defaults to \code{sqrt(.Machine$double.eps)}.}
    \item{\code{verbose}}{a logical indicating whether to provide
      some output on minimization progress.
      Defaults to \code{getOption("verbose")}.}
  }
  
  Additional control parameters for method \code{"genetic"} are as
  follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximum number of
      iterations for the genetic algorithm.  Defaults to 12.}
    \item{\code{popsize}}{an integer giving the population size for the
      genetic algorithm.  Default: 6.
      Only used if \code{start} is not given.}
    \item{\code{mutations}}{a number between 0 and 1 giving the
      probability of mutation per iteration.  Defaults to 0.1.}
  }

  Additional control parameters for method \code{"pclust"} are as
  follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      fixed-point iterations to be performed.  Default: 100.}
    \item{\code{nruns}}{an integer giving the number of fixed-point runs
      to be performed.  Default: 1.
      Only used if \code{start} is not given.}
    \item{\code{maxchains}}{an integer giving the maximal length of the
      Kernighan-Lin chains.  Default: 0 (no first variation improvements
      are attempted).}
  }
    
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

  Control parameters for method \code{"gmeans"} are as follows.
  \describe{
    \item{\code{gmeans}}{the path to the \code{gmeans} executable.}
    \item{\code{verbose}}{as for the genetic algorithm.}
    \item{\code{control}}{a character string specifying arguments passed
      over to the \code{gmeans} executable.}
  }

  Control parameters for method \code{"kmndirs"} are as follows.
  \describe{
    \item{\code{nstart}}{an integer giving the number of starting points
      to compute the starting value for the iteration stage.
      Default: 100.}
    \item{\code{maxiter}}{an integer giving the maximum number of
      iterations.
      Default: 10.}
  }

  Method \code{"CLUTO"} requires that the CLUTO \code{vcluster}
  executable is available.  CLUTO binaries for the Linux, Sun, OSX, and
  MS Windows platforms can be downloaded from
  \url{http://www-users.cs.umn.edu/~karypis/cluto/}.
  If the executable cannot be found in the system path via
  \code{\link{Sys.which}("vcluster")} (i.e., named differently or not
  made available in the system path), its (full) path must be specified
  in control option \code{vcluster}.

  Method \code{"gmeans"} requires that the \code{gmeans} executable is
  available.
  Sources for compilation with ANSI C++ compliant compilers
  are available from
  \url{http://www.dataminingresearch.com/index.php/2010/06/gmeans-clustering-software-compatible-with-gcc-4};
  original sources can be obtained from
  \url{http://userweb.cs.utexas.edu/users/dml/Software/gmeans.html};
  If the executable cannot be found in the system path via
  \code{\link{Sys.which}("gmeans")} (i.e., named differently or not
  made available in the system path), its (full) path must be specified
  in control option \code{vcluster}.
  
  Method \code{"kmndirs"} requires package \pkg{kmndirs} (currently, not
  available from CRAN), which provides an R interface to a suitable
  modification of the C code for the \eqn{k}-mean-directions algorithm
  made available as supplementary material to Maitra & Ramler (2010) at
  \url{http://pubs.amstat.org/doi/suppl/10.1198/jcgs.2009.08155}. 
  
  User-defined methods must have formals \code{x}, \code{k} and
  \code{control}, and optionally may have formals \code{weights} 
  or \code{m} if providing support for case weights or soft spherical
  \eqn{k}-means partitions, respectively.
}
\value{
  An object inheriting from classes \code{skmeans} and \code{pclust}
  (see the information on \link[clue:pclust]{pclust objects} in package
  \pkg{clue} for further details) representing the obtained spherical
  \eqn{k}-means partition, which is a list with components including the
  following:
  \item{prototypes}{a dense matrix with \code{k} rows giving the
    prototypes.}
  \item{membership}{cluster membership as a matrix with \code{k}
    columns (only provided if \eqn{m > 1}).}
  \item{cluster}{the class ids of the closest hard partition (the
    partition itself if \eqn{m = 1}).}
  \item{value}{the value of the criterion.}

  Objects representing spherical \eqn{k}-means partitions have special
  methods for \code{\link{print}},
  \code{\link[clue]{cl_validity}} (providing the \dQuote{dissimilarity
    accounted for}) from package \pkg{clue},
  and
  \code{\link[cluster]{silhouette}} from package \pkg{cluster} (the
  latter two take advantage of the special structure of the cosine
  distance to avoid computing full object-by-object distance matrices,
  and hence also perform well for large data sets).

  Package \pkg{clue} provides additional methods for objects inheriting
  from class \code{pclust}, see the examples.
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

  R. Maitra and I. P. Ramler (2010).
  A \eqn{k}-mean-directions algorithm for fast clustering of data on the
  sphere.
  \emph{Journal of Computational and Graphical Statistics}, \bold{19}/2,
  377--396.
}
\author{
  Kurt Hornik \email{Kurt.Hornik@wu.ac.at}, \cr
  Ingo Feinerer \email{feinerer@logic.at}, \cr
  Martin Kober \email{martin.kober@wu.ac.at}.
}
\examples{
set.seed(1234)

## Use CLUTO dataset 're0' and the reader for CLUTO sparse matrix
## format in package 'slam'.  (In text clustering applications, x will
## often be a DocumentTermMatrix object obtained from package 'tm'.)
x <- slam::read_stm_CLUTO(system.file("cluto", "re0.mat",
                                      package = "skmeans"))
## Which is not really small:
dim(x)

## Hard partition into 5 clusters.
hparty <- skmeans(x, 5, control = list(verbose = TRUE))
## Criterion value obtained:
hparty$value
## Compare with "true" classifications:
class_ids <- attr(x, "rclass")
table(class_ids, hparty$cluster)
## (Note that there are actually 10 "true" classes.)

## Plot the silhouette information for the obtained partition.
require("cluster")
plot(silhouette(hparty))
## Clearly, cluster 3 is "best", and cluster 5 needs splitting.

## Soft partition into 5 clusters.
sparty <- skmeans(x, 5, m = 1.1,
                  control = list(nruns = 5, verbose = TRUE))
## Criterion value obtained:
sparty$value
## (This should be a lower bound for the criterion value of the hard
## partition.)

## Compare the soft and hard partitions:
table(hparty$cluster, sparty$cluster)
## Or equivalently using the high-level accessors from package 'clue':
require("clue")
table(cl_class_ids(hparty), cl_class_ids(sparty))
## Which can also be used for computing agreement/dissimilarity measures
## between the obtained partitions.
cl_agreement(hparty, sparty, "Rand")

## How fuzzy is the obtained soft partition?
cl_fuzziness(sparty)
## And in fact, looking at the membership margins we see that the
## "sureness" of classification is rather high:
summary(cl_margin(sparty))
}
\keyword{cluster}
