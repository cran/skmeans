## Spherical (possibly sparse) k-means.
## Authors: Kurt Hornik, Ingo Feinerer, Martin Kober.

## Partition given vectors x_b by minimizing the criterion
##   \sum w_b u_{bj}^m d(x_b, p_j)
## where the w_b are case weights, u_{bj} is the membership of x_b to
## class j, p_j is the prototype of class j (and thus minimizes
##   \sum_b w_b u_{bj}^m d(x_b, p)
## over p), and d is *cosine dissimilarity*
##   d(x, p) = 1 - cos(x, p),
## where
##   cos(x, p) = <x, p> / (||x|| * ||p||)
## (so that d(x, p) is half the squared euclidean distance between the
## normalized x and the normalized p).
##
## (E.g., http://www.ams.jhu.edu/~priebe/.FILES/applepres.pdf or
## http://reference.wolfram.com/mathematica/ref/CosineDistance.html for
## similar use of cosine dissimilarity/distance.) 
##
## The literature only considers the unweighted hard spherical k-means
## problem (all w_b = 1 and m = 1), which we refer to as the "simple
## spherical k-means problem".  This has the criterion function
##   \sum_j \sum_{b in class j} d(x_b, p_j)
##     = B - \sum_j \sum_{b in class j} cos(x_b, p_j)
## with B the number of objects, so that minimizing the dissimilarity
## criterion is equivalent to maximizing
##   \sum_j \sum_{b in class j} cos(x_b, p_j))
## which is criterion I2 in CLUTO.
##
## All built-in simple spherical k-means methods internally use I2 but
## report the value of the corresponding dissimilarity criterion.

### * skmeans_family

## In the skmeans family object, we provide the "correct" cosine
## dissimilarity between objects and prototypes, irrespective of
## possible normalization.  For the methods, we internally always
## normalize objects and prototypes and hence use faster D and C
## functions.
## Prototypes are always represented by a dense matrix.
## Objects can be represented by dense matrices or sparse matrices
## provided that the following matrix computations work for these:
## * Subscripting rows
## * Scaling rows via multiplication (from the left) or division (from
##   the right) by a numeric vector
## * Taking (element-wise) squares
## * rowSums and colSums
## * tcrossprod with a dense matrix.

## Unfortunately, whereas the first three can easily be accomplished via
## (S3) subscript and Ops group methods, rowSums, colSums and tcrossprod
## are not generic, and S4 dispatch will not occur on the base function.
## Hence, we provide our own S3 generics and methods for supported
## classes.

g_tcrossprod <-
function(x, y = NULL)
    UseMethod("g_tcrossprod")
g_tcrossprod.default <-
function(x, y = NULL)    
    base::tcrossprod(x, y)
g_tcrossprod.simple_triplet_matrix <-
function(x, y = NULL)
    slam::tcrossprod_simple_triplet_matrix(x, y)
g_tcrossprod.dgCMatrix <-
function(x, y = NULL)
    Matrix::tcrossprod(x, y)
g_tcrossprod.dgTMatrix <-
function(x, y = NULL)
    Matrix::tcrossprod(x, y)

## (Using
##    g_tcrossprod.default <- base::tcrossprod
## has R CMD check NOTE .Internal calls, using
##    g_tcrossprod.dgTMatrix <- Matrix::tcrossprod
## loads Matrix dependencies.)

g_row_sums <-
function(x, na.rm = FALSE, dims = 1, ...) 
    UseMethod("g_row_sums")
g_row_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)     
    slam:::row_sums.simple_triplet_matrix(x, na.rm = na.rm,
                                         dims = dims, ...)
g_row_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::rowSums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)
g_row_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::rowSums(x, na.rm = na.rm, dims = dims, ...)

g_col_sums <-
function(x, na.rm = FALSE, dims = 1, ...)
    UseMethod("g_col_sums")
g_col_sums.simple_triplet_matrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    slam:::col_sums.simple_triplet_matrix(x, na.rm = na.rm,
                                         dims = dims, ...)
g_col_sums.default <-
function(x, na.rm = FALSE, dims = 1, ...)
    base::colSums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.dgCMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)
g_col_sums.dgTMatrix <-
function(x, na.rm = FALSE, dims = 1, ...)
    Matrix::colSums(x, na.rm = na.rm, dims = dims, ...)

skmeans_family <-
    clue::pclust_family(D =
                        function(x, prototypes) {
                            1 - g_tcrossprod(row_normalize(x),
                                             row_normalize(prototypes))
                        },
                        C =
                        function(x, weights, control) {
                            ## Computes weighted averages of the
                            ## normalized data.
                            x <- row_normalize(x)
                            weights <- weights / sum(weights)
                            out <- col_sums_with_weights(x, weights)
                            ## <NOTE>
                            ## Should prototypes be normalized?
                            ## They are not in the basic references.
                            ## If normalization is desired, uncomment
                            ##    out <- out / sqrt(sum(out ^ 2))
                            ## </NOTE>
                            out
                        },
                        init =
                        function(x, k) {
                            ## <FIXME>
                            ## Should perhaps ensure that this returns
                            ## unique prototypes.
                            as.matrix(x[sample.int(nrow(x), k), ,
                                        drop = FALSE])
                            ## </FIXME>
                        },
                        description = "spherical k-means")

### * skmeans

skmeans <-
function(x, k, method = NULL, m = 1, weights = 1, control = list())
{
    if(!all(row_norms(x) > 0))
        stop("Zero rows are not allowed.")
    
    ## Methods must at least have formals x, k and control.
    ## Try to ensure that formals m and weights are only used if
    ## different from the default ("simple") case.

    args <- list(x = x, k = k, control = control)
    if(!missing(m) && !identical(m, 1))
        args <- c(args, list(m = m))
    if(!missing(weights) && !all(weights == 1))
        args <- c(args,
                  list(weights = rep(weights, length.out = nrow(x))))

    skmeans_methods <-
        c(genetic = "genetic",
          pclust = "pclust",
          CLUTO = "CLUTO",
          LIH = "local_improvement_heuristic",
          LIHC = "local_improvement_heuristic_with_chains")

    if(!is.function(method)) {
        method <- if(is.null(method)) {
            if(is.null(args$m)) "genetic" else "plcust"
        } else if(is.character(method)) {
            pos <- pmatch(tolower(method),
                          tolower(names(skmeans_methods)))
            if(is.na(pos))
                stop(gettextf("Invalid skmeans method '%s'.", method))
            method <- skmeans_methods[pos]
        } else {
            stop("Invalid skmeans method.")
        }
        method <- get(sprintf(".skmeans_%s", method))
    }

    ## Now check whether the method has the required formals.
    na <- names(args)
    nf <- names(formals(method))
    if(any(ind <- is.na(match(na, nf))))
        stop(gettextf("Given skmeans method lacks formals %s",
                      paste(sQuote(na[ind]), collapse = " and ")))
    ## Otherwise, ensure that method gets called with "all" arguments so
    ## that it does not also have to provide defaults for these.
    if(("m" %in% nf) && !("m" %in% na))
        args <- c(args, list(m = m))
    if(("weights" %in% nf) && !("weights" %in% na))
        args <- c(args,
                  list(weights = rep(weights, length.out = nrow(x))))

    ## Call the skmeans method.
    do.call(method, args)
}
    
### * .skmeans_pclust

.skmeans_pclust <-
function(x, k, m, weights, control)
{
    ## Internally, we always keep x and prototypes normalized.
    .D_for_normalized_x_and_prototypes <-
        function(x, prototypes)
            1 - g_tcrossprod(x, prototypes)
    ## Add the 2 so that we do not have to recompute the value which is
    ## a nuisance if m > 1.
    .C_for_normalized_x_and_prototypes <-
        function(x, weights, control) {
            out <- col_sums_with_weights(x, weights)
            out / sqrt(sum(out ^ 2))
        }

    if(!is.null(start <- control$start)) {
        control$start <- if(inherits(start, "skmeans"))
            list(row_normalize(start$prototypes))
        else if(!is.list(start))
            list(row_normalize(start))
        else
            lapply(start, row_normalize)
    }

    x <- row_normalize(x)
            
    family <- skmeans_family
    C <- family$C
    D <- family$D
    family$C <- .C_for_normalized_x_and_prototypes
    family$D <- .D_for_normalized_x_and_prototypes
    out <- clue::pclust(x, k, family, m, weights, control)
    out$family$C <- C
    out$family$D <- D

    class(out) <- unique(c("skmeans", class(out)))
    out
}

### * .skmeans_genetic

.skmeans_genetic <-
function(x, k, weights = NULL, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 12L
    mutations <- control$mutations
    if(is.null(mutations))
        mutations <- 0.1
    popsize <- control$popsize
    if(is.null(popsize))
        popsize <- 6L
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-8
    space_for_time <- control$space_for_time
    if(is.null(space_for_time))
        space_for_time <- TRUE
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## In the weighted case, optimal prototypes are the weighted sums of
    ## the normalized objects in a class.  We can save time at the cost
    ## of space by precomputing the weighted normalized objects as wx.
    w <- weights
    if(all(w == 1)) {
        v <- nr
        wx <- x
        w <- NULL
    } else {
        v <- sum(w)
        if(space_for_time) {
            wx <- w * x
            w <- NULL
        } else {
            wx <- x
        }
    }

    ## Initialize p.
    p <- control$start
    if(inherits(p, "skmeans")) {
        p <- list(row_normalize(p$prototypes))
    } else if(!is.list(p)) {
        ## Initialize p by choosing k random rows of x.
        ## (Hence initial prototypes are already normalized.)
        p <- replicate(popsize,
                       as.matrix(x[sample.int(nr, k), , drop = FALSE]),
                       simplify = FALSE)
    }
    ## Initialize ids.
    ids <- lapply(p,
                  function(p1)
                  ids_from_similarities(g_tcrossprod(x, p1), k))
    
    if(verbose)
        message(gettextf("Initial solutions created (length: %d)",
                         length(ids)))
    
    value <- rep.int(0, popsize)

    genetic_mutator <- function(id) {
        ## Randomly change cluster membership .
        ## Simplistic solution.
        mutations2 <- mutations * k / (k - 1L) 
        newid <- id
        mut <- ceiling(runif(nr * mutations2, 0, nr))
        newid[mut] <- ceiling(runif(nr * mutations2, 0, k))
        newid
    }
    
    genetic_selection <- function(v) {
        v <- v + runif(length(v), min(v), max(v))
        which(v >= sort(v, decreasing = TRUE)[popsize])
    }

    iter <- 1L
    while(iter <= maxiter) {
        ## <FIXME>
        ## Make sure that we have no empty clusters.
        ## </FIXME>
        ## select best prototypes
        sel <- genetic_selection(value)
        if(verbose)
            message(gettextf("Selecting new population: %s",
                             paste(sel, collapse = ",")))
        ids <- ids[sel]
        p <- p[sel]
        value <- value[sel]
        for (i in seq_len(popsize)) {
            ## Generate new prototypes by mutating and k-means.
            new <- popsize + i
            ids[[new]] <- genetic_mutator(ids[[i]])
            sums <- .simple_skmeans_C_for_normalized_x(wx, ids[[new]],
                                                       k, nc, w)
            norms <- row_norms(sums)
            if(verbose)
                message(gettextf("Iteration: %d [KM] *** value[%d]: %g pre-optimized",
                                 iter, new, v - sum(norms)))
            repeat {
                p[[new]] <- sums / norms
                similarities <- g_tcrossprod(x, p[[new]])
                ids[[new]] <- ids_from_similarities(similarities, k)
                sums <- .simple_skmeans_C_for_normalized_x(wx, ids[[new]],
                                                           k, nc, w)
                norms <- row_norms(sums)
                oldvalue <- value[new]
                value[new] <- sum(norms)
                if(!is.na(oldvalue) && 
                    abs(oldvalue - value[new]) < reltol*(oldvalue + reltol)) break
            }
            if(verbose)
                message(gettextf("Iteration: %d [KM] *** value[%d]: %g",
                                 iter, new, v - value[new]))

        }
        iter <- iter + 1L
    }

    ## Get ids from the winner population.
    ids <- ids[[which.max(value)]]

    .simple_skmeans_object_for_normalized_x(wx, ids, k, nr, nc, w, v)
}

### * .skmeans_local_improvement_heuristic

.skmeans_local_improvement_heuristic <-
function(x, k, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-8
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## Initialize p.
    p <- .simple_skmeans_init_for_normalized_x(x, k, control)

    old_value <- 0

    iter <- 1L
    while(iter <= maxiter) {
        similarities <- g_tcrossprod(x, p)
        ids <- ids_from_similarities(similarities, k)
        ## New prototypes.
        sums <- .simple_skmeans_C_for_normalized_x(x, ids, k, nc)
        norms <- row_norms(sums)
        p <- sums / norms
        new_value <- sum(norms)
        if(verbose)
            message(gettextf("Iteration: %d [KM] *** value: %g",
                             iter, nr - new_value))
        if(abs(old_value - new_value)
           < reltol * (abs(old_value) + reltol)) {
            ## Block update performed only minor improvement.
            ## Try first variation steps.
            crossprods <- 2 * sweep(similarities, 2L, norms, "*")
            count <- 1L
            repeat {
                nids <- norms[ids]
                ind <- cbind(seq_len(nr), ids)
                ## The effects of removing an object from its cluster.
                new_norms_rem <- sqrt(nids ^ 2 - crossprods[ind] + 1)
                Delta_Q_rem <- new_norms_rem - nids
                ## The effects of adding an object to another cluster.
                new_norms_add <-
                    sqrt(sweep(crossprods, 2, norms ^ 2 + 1, "+"))
                Delta_Q_add <- sweep(new_norms_add, 2, norms, "-")
                ## What is the best such move?
                Delta_Q <- Delta_Q_rem + Delta_Q_add
                Delta_Q[ind] <- 0
                Delta_Q <- as.numeric(Delta_Q)
                pos <- which.max(Delta_Q)
                ## This changes object o from its cluster i (ids[o]) to
                ## cluster j.  If the change is large enough, perform
                ## it; otherwise we are done.
                if(Delta_Q[pos] < reltol * (abs(new_value) + reltol))
                    break
                j <- (pos - 1L) %/% nr + 1L
                o <- (pos - 1L) %% nr + 1L
                i <- ids[o]
                ids[o] <- j
                if(verbose)
                    message(gettextf("Moving %d from cluster %d to %d.",
                                     o, i, j))
                p[i, ] <- (sums[i, ] - c(x[o, ])) / new_norms_rem[o]
                p[j, ] <- (sums[j, ] + c(x[o, ])) / new_norms_add[o, j]
                new_value <- new_value + Delta_Q[pos]
                if(verbose)
                    message(gettextf("Iteration: %d [FV %d] *** value: %g",
                                     iter, count, nr - new_value))
                count <- count + 1L
                ## <FIXME>
                ## Of course, this is terribly inefficient, but let's
                ## just see if we can find better solutions by repeating
                ## first variation ...
                norms[c(i, j)] <-
                    c(new_norms_rem[o], new_norms_add[o, j])
                sums[c(i, j), ] <-
                    norms[c(i, j)] * p[c(i, j), ]
                crossprods[, c(i, j)] <-
                    2 * g_tcrossprod(x, sums[c(i, j), ])
                ## </FIXME>
            }
            if(count == 1L) break
        }
        old_value <- new_value
        iter <- iter + 1L
    }

    .simple_skmeans_object_for_normalized_x(x, ids, k, nr, nc)
}

### * .skmeans_local_improvement_heuristic_with_chains

.skmeans_local_improvement_heuristic_with_chains <-
function(x, k, control = NULL)
{
    maxiter <- control$maxiter
    if(is.null(maxiter))
        maxiter <- 100L
    maxchains <- control$maxchains
    if(is.null(maxchains))
        maxchains <- 10L
    reltol <- control$reltol
    if(is.null(reltol))
        reltol <- 1e-8
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")

    nr <- nrow(x)
    nc <- ncol(x)

    ## Normalize x.
    x <- row_normalize(x)

    ## Initialize p.
    p <- .simple_skmeans_init_for_normalized_x(x, k, control)

    old_value <- 0

    iter <- 1L
    while(iter <= maxiter) {
        similarities <- g_tcrossprod(x, p)
        ids <- ids_from_similarities(similarities, k)
        ## New prototypes.
        sums <- .simple_skmeans_C_for_normalized_x(x, ids, k, nc)
        norms <- row_norms(sums)
        p <- sums / norms
        new_value <- sum(norms)
        if(verbose)
            message(gettextf("Iteration: %d [KM] *** value: %g",
                             iter, nr - new_value))
        if(abs(old_value - new_value)
           < reltol * (abs(old_value) + reltol)) {
            ## Block update performed only minor improvement.
            ## Try first variation steps.
            crossprods <- 2 * sweep(similarities, 2L, norms, "*")
            count <- 1L
            repeat {
                ## Save prototypes for rollback.
                oldp <- p
                change <- 0L
                marked <- integer()
                chains <- 1L
                ## Try maxchains steps before discarding first variation
                ## (Kernighan-Lin chains)
                while (chains <= maxchains) {
                    nids <- norms[ids]
                    ind <- cbind(seq_len(nr), ids)
                    ## The effects of removing an object from its
                    ## cluster.
                    new_norms_rem <-
                        sqrt(nids ^ 2 - crossprods[ind] + 1)
                    Delta_Q_rem <- new_norms_rem - nids
                    ## The effects of adding an object to another
                    ## cluster.
                    new_norms_add <-
                        sqrt(sweep(crossprods, 2L, norms ^ 2 + 1, "+"))
                    Delta_Q_add <- sweep(new_norms_add, 2, norms, "-")
                    ## What is the best such move?
                    Delta_Q <- Delta_Q_rem + Delta_Q_add
                    Delta_Q[ind] <- -1e9
                    Delta_Q[marked,] <- -1e9
                    Delta_Q <- as.numeric(Delta_Q)
                    pos <- which.max(Delta_Q)
                    ## This changes object o from its cluster i (ids[o])
                    ## to cluster j.  If the change is large enough,
                    ## perform it; otherwise we are done.
                    j <- (pos - 1L) %/% nr + 1L
                    o <- (pos - 1L) %% nr + 1L
                    marked <- c(marked, o)
                    i <- ids[o]
                    ids[o] <- j
                    if(verbose)
                        message(gettextf("Moving %d from cluster %d to %d.",
                                         o, i, j))
                    p[i, ] <- (sums[i, ] - c(x[o, ])) / new_norms_rem[o]
                    p[j, ] <- (sums[j, ] + c(x[o, ])) / new_norms_add[o, j]
                    change <- Delta_Q[pos] + change
                    if(verbose)
                        message(gettextf("Iteration: %d [FV %d, CH %d] *** value: %g change: %g",
                                         iter, count, chains,
                                         nr - new_value, - 2 * change))
                    
                    ## <FIXME>
                    ## Of course, this is terribly inefficient, but
                    ## let's just see if we can find better solutions by
                    ## repeating first variation ...
                    norms[c(i, j)] <-
                        c(new_norms_rem[o], new_norms_add[o, j])
                    sums[c(i, j), ] <-
                        norms[c(i, j)] * p[c(i, j), ]
                    crossprods[, c(i, j)] <-
                        2 * g_tcrossprod(x, sums[c(i, j), ])
                    ## </FIXME>
                    if(change > reltol * (abs(new_value) + reltol)) {
                        count <- count + 1L
                        break	
                    }
                    chains <- chains + 1L	
                }
                if (chains == maxchains + 1L) {
                    p <- oldp
                    if (verbose) 
                        message(gettextf("Rolling back %d chain moves.",
                                         chains - 1L))
                    break
                }
                new_value <- sum(norms)
            }
            if(count == 1L) break
        }
        old_value <- new_value
        iter <- iter + 1L
    }

    ## Fix ids from chain evaluation.
    ids <- ids_from_similarities(g_tcrossprod(x, p), k)

    .simple_skmeans_object_for_normalized_x(x, ids, k, nr, nc)
}

### * .skmeans_CLUTO 

## Spherical sparse k-means via CLUTO file I/O.
    
.skmeans_CLUTO <-
function(x, k, control = NULL)
{
    ## Try to determine path to CLUTO vcluster executable.
    vcluster <- control$vcluster
    if(is.null(vcluster))
        vcluster <- "vcluster"
    if(Sys.which("vcluster") == "")
        stop("CLUTO vcluster executable not found")
    
    colmodel <- control$colmodel
    if(is.null(colmodel))
        colmodel <- "none"
    verbose <- control$verbose
    if(is.null(verbose))
        verbose <- getOption("verbose")
    ## Could add more fancy control list expansion eventually.
    control <- paste(as.character(control$control), collapse = " ")

    datfile <- tempfile()
    clustfile <- paste(datfile, ".clustering.", k, sep = "")
    on.exit(file.remove(datfile, clustfile))
    
    writeCM(x, datfile)
    cmdout <-
        system(sprintf("%s -colmodel=%s -clustfile=%s %s %s %d",
                       vcluster, colmodel, clustfile, control,
                       datfile, k),
               intern = TRUE)
    if(verbose)
        message(paste(cmdout, collapse = "\n"))
    
    result <- read.table(clustfile, header = FALSE) + 1
    ids <- result$V1

    .simple_skmeans_object_for_normalized_x(row_normalize(x), ids, k)
}

## * Helpers

## In the "simple" case (m = 1, identical weights) with normalized x:
## * the consensus function is the sum of all objects in a group;
## * the value is the sum of the norms of the consensus sums.

## Vectorized consensus for the simple normalized x case:

.simple_skmeans_C_for_normalized_x <- 
function(x, ids, k, nc = ncol(x), w = NULL)
{
    ## Only compute prototypes for used ids.
    all_ids_used <- sort(unique(ids))
    out <- matrix(0, length(all_ids_used), nc)
    if(!is.null(w)) {
        for(i in all_ids_used) {
            ## Save some space over using x <- w * x at once (which is
            ## the whole point of keeping w separate from x).
            out[i, ] <-
                g_col_sums(w[ids == i] * x[ids == i, , drop = FALSE])
        }
    } else {
        for(i in all_ids_used)
            out[i, ] <- col_sums_with_logical_index(x, ids == i)
    }
    rownames(out) <- all_ids_used
    out
}

## Initializer for the simple normalized x case.

.simple_skmeans_init_for_normalized_x <-
function(x, k, control)
{
    p <- control$start
    if(inherits(p, "skmeans")) {
        p <- p$prototypes
    } else {
        ## In case this is a list:
        if(is.list(p))
            p <- p[[1L]]
        ## In case this got us nowhere:
        if(is.null(p)) {
            ## Initialize p by choosing k random rows of x.
            ## (Hence initial prototypes are already normalized.)
            p <- as.matrix(x[sample.int(nrow(x), k), , drop = FALSE])
        }
    }
    p
}

## Object generator for the simple normalized x case.

.simple_skmeans_object_for_normalized_x <-
function(x, ids, k, nr = nrow(x), nc = ncol(x), w = NULL, v = NULL)
{
    sums <- .simple_skmeans_C_for_normalized_x(x, ids, k, nc, w)
    norms <- row_norms(sums)
    if(is.null(v))
        v <- if(is.null(w)) nr else sum(w)

    u <- clue::cl_membership(clue::as.cl_membership(ids), k)
    out <- clue::pclust_object(prototypes = sums / norms,
                               membership = u,
                               cluster = ids,
                               family = skmeans_family,
                               m = 1,
                               value = v - sum(norms))
    class(out) <- unique(c("skmeans", class(out)))
    out
}
    
## I2 for the simple case.

I2 <-
function(x, ids)
{
    x <- row_normalize(x)
    tab <- unique(ids)
    sums <- .simple_skmeans_C_for_normalized_x(x,
                                               match(ids, tab),
                                               length(tab),
                                               ncol(x))
    sum(row_norms(sums))
}

### * Utilities

row_norms <- 
function(x)
    sqrt(g_row_sums(x ^ 2))

row_normalize <-
function(x)
    x / row_norms(x)

col_sums_with_weights <-
function(x, w)
{
    ## Improve performance by leaving out Ops dispatch.
    ## (Other sparse matrix classes could be dealt with similarly ...)
    if(inherits(x, "simple_triplet_matrix")) {
        x$v <- x$v * w[x$i]
        slam:::col_sums.simple_triplet_matrix(x)
    } else if(inherits(x, "dgCMatrix")) {
        x@x <- x@x * w[x@i + 1L]
        Matrix::colSums(x)
    }
    else
        g_col_sums(w * x)
}

col_sums_with_logical_index <-
function(x, l)
{
    ## Improve performance by leaving out Ops dispatch.
    ## (Other sparse matrix classes could be dealt with similarly ...)
    if(inherits(x, "simple_triplet_matrix")) {
        x$v <- x$v * l[x$i]
        slam:::col_sums.simple_triplet_matrix(x)
    } else if(inherits(x, "dgCMatrix")) {
        x@x <- x@x * l[x@i + 1L]
        Matrix::colSums(x)
    }
    else
        g_col_sums(x[l, , drop = FALSE])
}

ids_from_similarities <-
function(x, k)
{
    ids <- max.col(x)
    all_ids_used <- sort(unique(ids))
    if(length(all_ids_used) < k) {
        ## Assign objects with the smallest similarities to their own
        ## cluster.
        o <- order(x[cbind(seq_along(ids), ids)])
        unused <- setdiff(seq_len(k), all_ids_used)
        ids[o[seq_along(unused)]] <- unused
    }
    ids
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
