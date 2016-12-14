
#' initializes data frame with discrete genes as rows and samples as columns.
#' each element is an integer representing how many qualifying variants fall
#' in that gene for the particular sample in question.
#' @param geneset set of genes to build collapsing matrix rows from.
#' @param samples set of samples to build collapsing matrix columns from.
#' @param qualvar qualvar data frame to code collapsing matrix on.
#' @export
CollapsingMatrixInit <- function(geneset, samples,
                                 qualvar=NULL) {
  mat.colnames <- samples
  mat.rownames <- geneset

  mat <- matrix(0, length(mat.rownames), length(mat.colnames))
  colnames(mat) <- mat.colnames
  rownames(mat) <- mat.rownames
  mat <- data.frame(mat, check.names=F)
  if (is.null(qualvar) == F) {
    mat <- CollapsingMatrixRecode(mat, qualvar)
  }
  return(mat)
}

#' read in a collapsing matrix file from file (gzip detection supported) into
#' CollapsingMatrix data frame.
#' @param mat.filename name of collapsing matrix file.
#' @export
CollapsingMatrixRead <- function(mat.filename) {
  mat <- ReadLargeTable(mat.filename, header=T, sep="\t", check.names=F)
  rownames(mat) <- mat[,1]
  mat[,1] <- NULL
  return(mat)
}

#' take a subset of collapsing matrix where each gene has at least one 
#' sample with a qualifying variant and each sample has at least one 
#' qualifying variant in a gene.
#' @param mat collapsing matrix.
#' @export
CollapsingMatrixShrink <- function(mat) {
  mat <- mat[rowSums(mat) > 0, , drop=F]
  samples <- colnames(mat)
  for (sample in samples) {
    if (max(mat[[sample]]) == 0) {
      mat[[sample]] <- NULL
    }
  }
  return(mat)
}

#' recode a collapsing matrix using a qualvar data frame.
#' @param mat collapsing matrix data frame.
#' @param qualvar qualvar data frame.
#' @export
CollapsingMatrixRecode <- function(mat, qualvar) {
  genes <- sort(unique(qualvar[,1]))
  # only use subset of qualvar where samples are in mat
  qualvar <- qualvar[ qualvar[,3] %in% colnames(mat), , drop=F]
  for (gene in genes) {
    qualvar.g <- qualvar[qualvar[,1] == gene, ,drop=F]
    sampleIDs.qualvar.g <- unique(qualvar.g[,3])
    for (sampleID in sampleIDs.qualvar.g) {
      qualvar.g.s <- qualvar.g[ qualvar.g[,3] == sampleID, , drop=F]
      mat[gene, sampleID] <- sum(qualvar.g.s[,4])
    }
  }
  return(mat)
}

#' write collapsing matrix to output filename (gzipping supported).
#' @param mat collapsing matrix data frame.
#' @param filename output filename.
#' @export
CollapsingMatrixWrite <- function(mat, filename) {
  mat.out <- data.frame("sample/gene"=rownames(mat),
                        mat, check.names=F)
  WriteLargeTable(mat.out, filename,
                  col.names=T, row.names=F,
                  quote=F, sep="\t")
}
