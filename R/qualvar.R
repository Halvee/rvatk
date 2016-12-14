#' read a file of format qualvar into a qualvar dataframe.
#' @param qualvar.filename name of qualvar file.
#' @param samples vector of samples to subset qualvar on.
#' @param geneset vector of genes to subset qualvar on.
#' @export
QualvarRead <- function(qualvar.filename,
                        samples = NULL,
                        geneset = NULL) {
  qv <- ReadLargeTable(qualvar.filename, header=F)
  qv <- QualvarSubset(qv, samples, geneset)
}

#' subset qualvar dataframe on user-defined sets of samples and genes.
#' @param qv qualvar data frame.
#' @param samples vector of samples to subset qualvar on.
#' @param geneset vector of genes to subset qualvar on.
QualvarSubset <- function(qv,
                          samples = NULL,
                          geneset = NULL) {
  if (is.null(samples) == F) {
    qv <- qv[ qv[,3] %in% samples, , drop=F]
  }
  if (is.null(geneset) == F) {
    qv <- qv[ qv[,1] %in% geneset, , drop=F]
  }
  return(qv)
}

#' convert qualvar dataframe to collapsing matrix data frame.
#' @param qv qualvar data frame.
#' @param samples vector of samples to subset qualvar on.
#' @param geneset vector of genes to subset qualvar on.
QualvarToCollapsingMatrix <- function(qv,
                                      samples,
                                      geneset) {
  qv <- QualvarSubset(qv, samples, geneset)
  clps <- CollapsingMatrixInit(geneset, samples, qualvar=qv)
  return(clps)
}

