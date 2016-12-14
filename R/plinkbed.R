#' functions for working with plinkbed, a snpMatrix object courtesy of the
#' snpStats pkg.
#' 
#' initialize plinkbed.
#' @param bed.file plink bed filename.
#' @param ... variables to be passed to the read.plink function in snpStats.
#' @importFrom snpStats read.plink
#' @export
PlinkBedInit <- function(bed.file,
                         ...) {
  pb <- read.plink(bed.file, ...)
  return(pb)
}

#' subset plinkbed on user-defined set of samples and variants.
#' @param pb plinkbed, a snpMatrix object courtesy of the snpStats pkg
#' @param samples vector of samples to subset plinkbed on.
#' @param vars vector of variants to subset plinkbed on.
#' @export
PlinkBedSubset <- function(pb,
                           samples = NULL,
                           vars = NULL) {
  if (is.null(samples) == F) {
    pb$genotypes <- pb$genotypes[samples, ]
    pb$fam <- pb$fam[samples, ]
  }
  if (is.null(vars) == F) {
    pb$genotypes <- pb$genotypes[, vars]
    pb$map <- pb$map[vars, ]
  }
  return(pb)
}

#' get sum of allele counts across plinkbed for user-defined 
#' samples and variants.
#' @param pb plinkbed, a snpMatrix object courtesy of the snpStats pkg
#' @param samples vector of samples to subset on.
#' @param vars vector of variants to subset on.
PlinkBedAlleleSum <- function(pb,
                              samples=NULL,
                              vars=NULL) {
  pb.subset <- PlinkBedSubset(pb,
                              samples,
                              vars)
  gts.transform <- PlinkBedGtsTransform(pb.subset$genotypes)
  gts.sum <- sum( gts.transform )
  return(gts.sum)
}

#' return variant calls per sample that match user-defined genotype.
#' @param pb plinkbed, a snpMatrix object courtesy of the snpStats pkg
#' @param samples a vector of samples to subset on.
#' @param vars a vector of variants to subset on.
#' @param gt genotype count to subset on.
#' @param as.df return calls in data.frame format?
#' @importFrom methods as
#' @export
PlinkBedCalls <- function(pb,
                          samples,
                          vars,
                          gt=1,
                          as.df = F) {
  gt.trans <- PlinkBedGtsTransform(gt)
  pb.subset <- PlinkBedSubset(pb,
                              samples,
                              vars)
  gts <- as(pb.subset$genotypes, "numeric")
  if (is.na(gt)) {
    calls.ij <- which(is.na(gts),
                      arr.ind = T,
                      useNames=F)
  } else {
    calls.ij <- which(gts==gt.trans,
                      arr.ind = T,
                      useNames=F)
  }
  calls.ij[,1] <- samples[ as.numeric(calls.ij[,1]) ]
  calls.ij[,2] <- vars[ as.numeric(calls.ij[,2]) ]
  if (as.df == T) {
    return(data.frame(calls.ij[,c(2,1),drop=F], stringsAsFactors=F))
  } else {
    return(calls.ij[,c(2,1),drop=F])
  }
}

#' return genotypes from plinkbed into counts that are expected.
#' @param gts vector of genotypes from plinkbed snpMatrix obj.
#' @export
PlinkBedGtsTransform <- function(gts) {
  return( -( as(gts, "numeric") - 2) )
}

#' convert plinkbed data to qualvar data frame and return.
#' @param pb plinkbed, a snpMatrix obj courtesy of the snpStats pkg.
#' @param samples a vector of samples to subset on.
#' @param varinfo a data frame in the varinfo format (varID,geneName,varClass).
#' @param varinfo.class variant class to subset on in varinfo data frame.
#' @export
PlinkBedToQualvar <- function(pb,
                              samples,
                              varinfo,
                              include.na = TRUE,
                              varinfo.class=NULL) {

  varinfo <- varinfo[,1:2]
  variants <- unique(varinfo[,1])
  pb.calls.ij.ac1 <- PlinkBedCalls(pb,
                                   samples,
                                   variants,
                                   gt=1,
                                   as.df=T)
  if (nrow(pb.calls.ij.ac1) > 0) {
    pb.calls.ij.ac1[,3] <- 1
    pb.calls.ij.ac1.g <- merge(pb.calls.ij.ac1,
                               varinfo, by=1)
  } else {
    pb.calls.ij.ac1.g <- data.frame(matrix(NA,nrow=0,ncol=4))
  }
  pb.calls.ij.ac2 <- PlinkBedCalls(pb,
                                   samples,
                                   variants,
                                   gt=2,
                                   as.df=T)
  if (nrow(pb.calls.ij.ac2) > 0) {
    pb.calls.ij.ac2[,3] <- 2
    pb.calls.ij.ac2.g <- merge(pb.calls.ij.ac2,
                               varinfo, by=1)
  } else {
    pb.calls.ij.ac2.g <- data.frame(matrix(NA,nrow=0,ncol=4))
  }
  if (include.na == T) {
    pb.calls.ij.acNA <- PlinkBedCalls(pb,         
                                      samples,                                 
                                      variants,                             
                                      gt=NA,                                    
                                      as.df=T)                                 
    if (nrow(pb.calls.ij.acNA) > 0) {                           
      pb.calls.ij.acNA[,3] <- NA 
      pb.calls.ij.acNA.g <- merge(pb.calls.ij.acNA,         
                                  varinfo, 
                                  by=1)
    } else {
      pb.calls.ij.acNA.g <- data.frame(matrix(NA,nrow=0,ncol=4))
    }
    pb.calls.ij.g <- rbind(pb.calls.ij.acNA.g, 
                           pb.calls.ij.ac1.g, 
                           pb.calls.ij.ac2.g) 
  } else {
    pb.calls.ij.g <- rbind(pb.calls.ij.ac1.g, pb.calls.ij.ac2.g)
  }
  pb.calls.ij.g <- pb.calls.ij.g[,c(4,1,2,3)]
  return(pb.calls.ij.g)
}
