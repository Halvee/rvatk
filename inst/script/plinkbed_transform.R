#!/usr/bin/Rscript

require(optparse)
require(data.table)

require(rvatk)                                                               

# GLOBAL VARS
OUTPUT.TYPES <- c("qualvar","long","collapsing.matrix")

Main <- function() {
  option.list <- GetOptList()
  opt <- parse_args(OptionParser(option_list=option.list))

  # read in varinfo file
  varinfo <- read.table(opt$varinfo.file, stringsAsFactors=F)
  if (ncol(varinfo) > 2 & is.null(opt$varinfo.class) == F) {
    varinfo <- varinfo[ varinfo[,3] == opt$varinfo.class, , drop=F]
  } else {
    varinfo <- unique(varinfo[,c(1,2)])
  }
  variants <- unique(varinfo[,1])

  # read in list of samples
  sampleped <- read.table(opt$sampleped.file, stringsAsFactors=F, sep="\t")
  cases <- sampleped[sampleped[,6]==2, 2]
  ctrls <- sampleped[sampleped[,6]==1, 2]
  samples <- c(cases, ctrls)
  
  # verify no repeats of sampleID
  if (TRUE %in% duplicated(samples)) {
    stop("duplicated samples in input sample file")
  }
  
  # initialize PlinkBed
  pb <- PlinkBedInit(opt$plink.bed)

  # subset varinfo on vars actually in PlinkBed
  varinfo <- varinfo[ varinfo[,1] %in% rownames(pb$map), ]

  qualvar <- PlinkBedToQualvar(pb,
                               samples,
                               varinfo, 
                               varinfo.class=opt$varinfo.class)
  qualvar <- qualvar[order(qualvar[,1], qualvar[,2], qualvar[,3]), ]
  # qualvar <- gene,var,sample,gt
  if (opt$output.type == "collapsing.matrix") {
    if (is.null(opt$geneset.file)) {
      cat("to create collapsing matrix output need geneset.\n")
      q()
    }
    # read in geneset file
    geneset <- ReadListFile(opt$geneset.file)

    # verify no duplicated genes in geneset
    if (TRUE %in% duplicated(geneset)) {
      stop("duplicated gene names in input geneset file")
    }

    # build collapsing matrix based on geneset and sample files
    mat <- CollapsingMatrixInit(geneset, samples)
    mat <- CollapsingMatrixRecode(mat, qualvar)                                    
    CollapsingMatrixWrite(mat, opt$output.file)
  } else if (opt$output.type == "long") {
    output <- qualvar
  } else {
    WriteLargeTable(qualvar, opt$output.file,
                    col.names=F,
                    row.names=F, quote=F, sep="\t")
  }
}

GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--plink-bed"), action="store", type="character",
                default=NULL, dest="plink.bed",
                help="plink .bed file,,
                      required, [default %default]"),

    make_option(c("--varinfo-file"), action="store", type="character", 
                default=NULL, dest="varinfo.file",
                help="varinfo file, required, [default %default]"),

    make_option(c("--varinfo-class"), action="store", type="character",
                default=NULL, dest="varinfo.class",
                help="variant class in varinfo file to subset on,
                      [default %default]"),

    make_option(c("--output-file"),
                action="store", type="character",
                default=NULL,
                dest="output.file",
                help="output file name, required, [default %default]"),

    make_option(c("--output-type"),
                action="store", type="character",
                default="collapsing.matrix",
                dest="output.type",
                help="format for output file, [default %default]"),

    make_option(c("--sampleped-file"), action="store", type="character",
                default=NULL, dest="sampleped.file", 
                help="name of required sample file, [default %default]"),

    make_option(c("--geneset-file"), action="store", 
                type="character",      
                default=NULL, 
                dest="geneset.file",                          
                help="name of required geneset file, [default %default]"), 

    make_option(c("--condition-file"), action="store", type="character",
                default=NULL, dest="condition.file",
                help="name of condition file to subset genotype csv on,
                      if no qualvar file is provided [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
  Main()
}
