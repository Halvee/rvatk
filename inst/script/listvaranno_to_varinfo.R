#!/usr/bin/Rscript

require(optparse)
require(data.table)

require(devtools)
devtools::load_all("/Users/mh3534/devel/rvatk")
#require(rvatk)

Main <- function() {
  option.list <- GetOptList()
  optarg <- parse_args(OptionParser(option_list=option.list),
                        positional_arguments = TRUE)
  opt <- optarg$options
  cnd.filenames <- optarg$args
  
  # read in comma-delim list of classes
  classes <- StringToSet(opt$condition.classes)

  # only continue if number of condition files matches number of class names
  if((length(classes) != length(cnd.filenames)) | 
     (length(classes) == 0) | (length(cnd.filenames) == 0)) {
    stop("number of class names not equal to number of condition files:
          run script as follows:
          listvaranno_to_varinfo <--condition-classes A,B,C> [OTHER_OPTS]
          <condition_file_A> <condition_file_B> <condition_file_C>")
  }

  # read in listvaranno csv file                                                
  lva <- ReadLargeTable(opt$listvaranno.csv, 
                        showProgress=F, 
                        header=T, check.names=F)
  print(lva[[kGeneNameCol]])
  lva[[kGeneNameCol]] <- gsub("'","",lva[[kGeneNameCol]])   

  var.names <- c()                                                              
  var.genes <- c()                                                              
  var.classes <- c() 
  for (i in 1:length(classes)) {
    cls <- classes[i]
    cnd.filename <- cnd.filenames[i]
    
    # read in condition file into table
    cnd <- ReadConditionFile(cnd.filename)

    # exclude 'loo maf' from cnd file if it is present, can't use here
    if (opt$include.loo.maf == FALSE) {
      cnd <- cnd[cnd[,1] != "loo maf", , drop=F]
    }

    # subset lva table on conditions defined in condition file
    lva.cls <- TblCndSubset(lva, cnd)
    
    # keep track of each variant, the gene it maps to and its variant class,
    # as defined by the condition file
    var.names <- c(var.names, lva.cls[[kVariantIdCol]])
    var.genes <- c(var.genes, lva.cls[[kGeneNameCol]])
    
    var.classes <- c(var.classes, rep(cls, nrow(lva.cls)) )
  }

  out.df <- data.frame(var.name = var.names, 
                       var.gene = var.genes,
                       var.class = var.classes)
  out.df <- out.df[order(out.df$var.name), ]

  write.table(out.df, file = opt$varinfo.file,
              row.names=F, col.names=F, quote=F, sep="\t")
}


GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--listvaranno-csv"), action="store", type="character",
                default=NULL, dest="listvaranno.csv", 
                help="name of ATAV listvaranno csv, 
                      [default %default], required"),

    make_option(c("--condition-classes"), action="store", type="character",
                default=NULL, dest="condition.classes",
                help="comma-delimited condition classes to subset variants on,
                      [default %default], required"),

    make_option(c("--varinfo-file"), action="store", type="character",
                default=NULL, dest="varinfo.file",
                help="name of output table listing gene and 
                      function class for each variant,
                      [default %default], required"),

    make_option(c("--include-loo-maf"), action="store_true", default=FALSE,      
                    dest="include.loo.maf",
                    help="do not remove 'loo maf' from condition table,
                          if it is present")
  )

  return(option.list)
}


if (interactive() == F) {
  Main()
}
