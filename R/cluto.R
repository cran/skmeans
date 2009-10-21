## Authors: Ingo Feinerer, Kurt Hornik

## Read in CLUTO (sparse) matrix format

readCM <-
function(con, clabel = NULL)
{
    ## Read in the matrix file
    l <- strsplit(readLines(con, warn = FALSE), "[[:space:]]+")
    l <- lapply(l, as.double)
    l <- lapply(l, na.omit)

    ## Extract header information
    nRow <- as.integer(l[[1L]][1L])
    nCol <- as.integer(l[[1L]][2L])
    nElem <- l[[1L]][3L]
    ## Remove header
    l <- l[-1L]

    ## Compute i, j, and v slots for a simple_triplet_matrix
    rowLen <- sapply(l, length)
    l <- unlist(l)
    i <- rep.int(seq_len(nRow), rowLen / 2)
    j <- l[seq.int(1, length(l), by = 2)]
    v <- l[seq.int(2, length(l), by = 2)]

    ## Sanity check
    if(nElem != length(v))
        stop("invalid matrix format")

    ## Generate sparse matrix
    mat <- simple_triplet_matrix(i, j, v, nRow, nCol)

    ## Use column labels file if available
    if(!is.null(clabel))
        colnames(mat) <- readLines(clabel)

    mat
}

## Write CLUTO (sparse) matrix format

writeCM <-
function(mat, con)
{
    mat <- as.simple_triplet_matrix(mat)

    ## Generate header
    header <- paste(nrow(mat), ncol(mat), length(mat$v))

    ## Generate content
    content <- Map(function(u, v) paste(u, v, collapse = " "),
                   split(mat$j, mat$i),
                   split(mat$v, mat$i))

    ## Write out
    writeLines(c(header, unlist(content)), con)
}
