#!/usr/bin/Rscript

require(optparse)
require(data.table)

require(rvatk)                                                               

# GLOBAL VARS
kOutputTypes <- c("qualvar","long","collapsing.matrix")

Main <- function() {
  option.list <- GetOptList()
  opt <- parse_args(OptionParser(option_list=option.list))

  # read in varinfo file
  varinfo <- read.table(opt$varinfo.file, stringsAsFactors=F)
  # if geneset passed, subset varinfo on geneset
  if (is.null(opt$geneset.file) == F) {
    geneset <- ReadListFile(opt$geneset.file)
    varinfo <- varinfo[ varinfo[,2] %in% geneset, , drop=F]
  }

  # read in list of samples
  sampleped <- read.table(opt$sampleped.file, stringsAsFactors=F, sep="\t")
  cases <- sampleped[sampleped[,6]==2, 2]
  ctrls <- sampleped[sampleped[,6]==1, 2]
  samples <- c(cases, ctrls)
  iid.fid <- SamplepedIidFidList(sampleped)
  
  # verify no repeats of sampleID
  if (TRUE %in% duplicated(samples)) {
    stop("duplicated samples in input sample file")
  }

  if (is.null(opt$plink.bed) == F) {
    # initialize PlinkBed
    pb <- PlinkBedInit(opt$plink.bed)
    
    # only keep varinfo vars that are in plinkbed
    varinfo <- varinfo[ varinfo[,1] %in% rownames(pb$map), ]

    qualvar.full <- PlinkBedToQualvar(pb,                                            
                                      samples,                                       
                                      varinfo)

  } else if (is.null(opt$qualvar) == F) {
    qualvar.full <- QualvarRead(opt$qualvar)
  } else {
    cat("ERROR: need to pass plink .bed or .qualvar file\n")
    q()
  }

  #
  var.classes <- sort(unique(varinfo[,3]))

  # build ouptut table
  tbl <- data.frame(matrix(0, 
                           nrow=length(unique(varinfo[,2])),
                           ncol=(1+2*length(var.classes))
                          )
                   )

  genes <- sort(unique(varinfo[,2]))
  var.classes <- sort(unique(varinfo[,3]))

  tbl <- data.frame(matrix(0,
                           nrow = length(genes),
                           ncol = length(var.classes)*2 + 1) )
  colnames(tbl)[1] <- "gene"
  tbl[,1] <- genes
  i <- 2
  for (var.class in var.classes) {
    colnames(tbl)[i] <- paste(var.class,"trans",sep=".")
    colnames(tbl)[i+1] <- paste(var.class,"ntrans",sep=".")
    i <- i + 2
  }
  rownames(tbl) <- genes
  
  for (var.class in var.classes) {
    varinfo.c <- varinfo[ varinfo[,3] == var.class, ]
    qualvar <- qualvar.full[ qualvar.full[,2] %in% varinfo.c[,1], , drop=F]
    qualvar <- qualvar[ qualvar[,4] == 1, ]
    for (gene in genes) {
      trans <- c()
      ntrans <- c()
      var.na <- c()
      qualvar.c.ctrl <- qualvar[ qualvar[,1] == gene &
                                 qualvar[,3] %in% ctrls, , drop=F]
      var.counts <- table(qualvar.c.ctrl[,2])
      if (nrow(qualvar.c.ctrl) > 0) {
        for (i in 1:nrow(qualvar.c.ctrl)) {
          variant <- qualvar.c.ctrl[i, 2]
          ctrlID <- qualvar.c.ctrl[i, 3]
          if (var.counts[[variant]] >= (opt$max.maf * length(ctrls))) { 
            next 
          }
          
          fid <- iid.fid[[ctrlID]]
          fid.cases <- sampleped[ sampleped[,1] == fid & sampleped[,6]==2, 2]
          for (fid.case in fid.cases) {
            qualvar.c.case <- qualvar[ which(qualvar[,2] == variant &
                                             qualvar[,3]  == fid.case), 
                                             , 
                                             drop = F ]
            if (nrow(qualvar.c.case) > 0) {
              if (is.na(qualvar.c.case[1, 4]) == T) {
                var.na <- c(var.na, fid.case)
              } else {
                trans <- c(trans, fid.case)
              }
            } else {
              
              ntrans <- c(ntrans, fid.case)
            }
          }
        } 
      }
      tbl[gene, paste(var.class, "trans", sep=".")] <- length(unique(trans))
      tbl[gene, paste(var.class, "ntrans", sep=".")] <- length(unique(ntrans))
    }
  }

  tbl <- tbl[rowSums(tbl[,-(1)]) > 0, ]

  write.table(tbl, file=opt$output.file, row.names=F, col.names=T,
              quote=F, sep="\t")
}

GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--plink-bed"), action="store", type="character",
                default=NULL, dest="plink.bed",
                help="plink .bed file, [default %default]"),

    make_option(c("--qualvar"), action="store", type="character", 
                default=NULL, dest="qualvar",
                help="qualvar file, [default %default]"),

    make_option(c("--varinfo-file"), action="store", type="character", 
                default=NULL, dest="varinfo.file",
                help="varinfo file, required, [default %default]"),

    make_option(c("--output-file"),
                action="store", type="character",
                default=NULL,
                dest="output.file",
                help="output file name, required, [default %default]"),

    make_option(c("--sampleped-file"), action="store", type="character",
                default=NULL, dest="sampleped.file", 
                help="name of required sample file, [default %default]"),

    make_option(c("--max-maf"), action="store", type="double",
                default=0.01, dest="max.maf",
                help="maximum maf allowed in unaffecteds, [default %default]"),

    make_option(c("--geneset-file"), action="store", 
                type="character",      
                default=NULL, 
                dest="geneset.file",                          
                help="name of required geneset file, [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
  Main()
}
