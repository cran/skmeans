## Write Gmeans (sparse) matrix format

# Gmeans uses a compressed column storage (CCS)
# See http://www.cs.utexas.edu/users/inderjit/Resources/sparse_matrices
#
# However since Gmeans clusters along columns, and the input for our skmeans clusters along rows,
# we would need to transpose the matrix first, and then write it to CCS.
#
# Instead we could directly write to compressed row storage (CRS) to avoid the transpose
# See http://netlib.org/linalg/html_templates/node92.html#SECTION00931200000000000000
writeGM <-
function(mat, file)
{
    x <- t(as.simple_triplet_matrix(mat))
    # Based on slam/work/Matrix.R
    ind <- order(x$j, x$i)
    write(paste(nrow(x), ncol(x), length(x$v)),
          sprintf("%s_dim", file))
    write(x$i[ind] - 1L,
          sprintf("%s_row_ccs", file), sep = "\n")
    write(c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
          sprintf("%s_col_ccs", file), sep = "\n")
    write(x$v[ind],
          sprintf("%s_tfn_nz", file), sep = "\n")
}
