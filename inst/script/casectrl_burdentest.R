#!/usr/bin/Rscript

require(optparse)
require(data.table)

require(rvatk)                                                               

kCasectrlAllColnames <- c("name","method","alternative")
kCasectrlFetColnames <- c(kCasectrlAllColnames,
                          "n.case.qual", "n.case.nonqual",
                          "n.ctrl.qual", "n.ctrl.nonqual",
                          "odds.ratio",
                          "conf.int.lower", "conf.int.upper",
                          "p.value")
kCasectrlPoissonColnames<- c(kCasectrlAllColnames,
                             "n.var.case", "n.case", "var.rate.case",
                             "n.var.ctrl", "n.ctrl", "var.rate.ctrl",
                             "conf.int.lower", "conf.int.upper",
                             "p.value")

Main <- function() {
  option.list <- GetOptList()
  opt <- parse_args(OptionParser(option_list=option.list))

  # read in sample file
  samples.tbl <- read.table(opt$sampleped.file, stringsAsFactors=F)
  cases <- samples.tbl[samples.tbl[,6]==2, 2]
  ctrls <- samples.tbl[samples.tbl[,6]==1, 2]
  samples <- c(cases, ctrls)

  # verify no repeats of sampleID or genename
  if (TRUE %in% duplicated(samples)) {
    stop("duplicated samples in input sample file")
  } 

  # read in qualvar file if defined and convert to matrix, 
  # else read in collapsing matrix
  if (is.null(opt$collapsing.matrix) == F) {
    mat <- CollapsingMatrixRead(opt$collapsing.matrix)
    
    # exclude genes where no cases or controls have qualifying variants           
    #mat <- mat[rowSums(mat) > 0, , drop=F]                                        
    mat <- CollapsingMatrixShrink(mat)

    geneset <- rownames(mat)
    if (length(geneset) > 0) {
      mat <- mat[geneset, , drop=F]
    }
  } else if (is.null(opt$qualvar) == F) {
    geneset <- ReadListFile(opt$geneset.file)
    qualvar <- QualvarRead(opt$qualvar,
                           samples = samples,
                           geneset = geneset)
    mat <- CollapsingMatrixInit(geneset, samples) 
    mat <- CollapsingMatrixRecode(mat, qualvar) 
    print(mat)
  } else {
    cat("require a qualvar file or a collapsing matrix file. Exiting...\n")
    q()
  }

  genesets <- c()
  tbl.nrow <- length(geneset)
  if (is.null(opt$geneset.file) == F) { 
    geneset <- scan(opt$geneset.file, what=character(), quiet=T) 
    geneset <- intersect(geneset, rownames(mat))  
    tbl.nrow <- length(geneset)
  } else if (is.null(opt$genesets.file) == F) {
    genesets <- ReadGenesetsFile(opt$genesets.file)
    for (geneset.name in names(genesets)) {
      genesets[[geneset.name]] <- intersect(genesets[[geneset.name]],
                                            rownames(mat) )
    }
    tbl.nrow <- length(names(genesets))
  }



  if (opt$collapse.thresh > 0) {
    res.tbl.colnames <- kCasectrlFetColnames
  } else {
    res.tbl.colnames <- kCasectrlPoissonColnames
  }
  res.tbl <- data.frame(matrix(0, nrow=tbl.nrow,
                               ncol=length(res.tbl.colnames)))
  colnames(res.tbl) <- res.tbl.colnames

  if (length(genesets) > 0) {
    geneset.names <- sort(names(genesets))
    for (i in 1:length(geneset.names)) {
      geneset <- genesets[[ geneset.names[i] ]]
      res <- BurdenTest(mat, geneset, cases, ctrls,
                        collapse.thresh=opt$collapse.thresh,
                        alternative=opt$alternative) 
      res$name <- geneset.names[i]
      res.tbl <- UpdateTable(res.tbl, i, res)       
    }
  }

  else if (opt$combine.geneset.stats == T) {
    res <- BurdenTest(mat, geneset, cases, ctrls,
                      collapse.thresh=opt$collapse.thresh,
                      alternative=opt$alternative)
    res$name <- "ALL"
    res.tbl <- UpdateTable(res.tbl, 1, res)
    res.tbl <- res.tbl[1, , drop=F]
  } else {
    for (i in 1:length(geneset)) {
      gene <- geneset[i]
      print(gene)
      # get gene/case, gene/ctrl subsets of matrix
      res <- BurdenTest(mat, gene, 
                        cases, ctrls,
                        collapse.thresh=opt$collapse.thresh,
                        alternative=opt$alternative)
      res$name <- gene
      res.tbl <- UpdateTable(res.tbl, i, res)
    }
  }
  write.csv(res.tbl, file=opt$out.csv, row.names=F, quote=F)
}


GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--out-csv"), action="store", type="character",
                default="out.csv", dest="out.csv",
                help="name of output csv, [default %default]"),

    make_option(c("--sampleped-file"), action="store", type="character",
                default=NULL, dest="sampleped.file", 
                help="name of required sample file, [default %default]"),

    make_option(c("--geneset-file"), action="store", type="character",      
                  default=NULL, dest="geneset.file",                          
                  help="name of required geneset file, [default %default]"), 

    make_option(c("--genesets-file"), action="store", type="character", 
                default=NULL, dest="genesets.file",
                help="name of required genesets file, [default %default]"),

    make_option(c("--qualvar"), action="store", type="character",
                default=NULL, dest="qualvar",
                help="name of .qualvar file, [default %default]"),

    make_option(c("--collapsing-matrix"), action="store", type="character",
                default=NULL, dest="collapsing.matrix",
                help="name of collapsing matrix file, 
                      if no qualvar file provided [default %default]"),

    make_option(c("--collapse-thresh"), action="store", type="integer",
                default=0, dest="collapse.thresh",
                help="minimum number of vars in gene for sample to qual,
                      set to 0 to instead use rate-based test,
                      [default %default]"),

    make_option(c("--combine-geneset-stats"), action="store_true", default=F,
                dest="combine.geneset.stats", 
                help="perform single test on all stats in 
                      user-defined geneset, [default %default]"),

    make_option(c("--alternative"), action="store", type="character", 
                default="two.sided", dest="alternative",
                help="alternative hypothesis for case/control test, [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
  Main()
}
