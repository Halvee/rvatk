#!/usr/bin/Rscript

### Libraries -----------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

### Functions -----------------------------------------------------------------

main <- function() {
  # Main loop of script.
  #
  # Returns:
  #   None.
  
  args <- commandArgs(trailing=T)

  option.list <- GetOptList()
  opt <- parse_args(OptionParser(option_list=option.list),
                    args=args)
  if (is.null(opt$listvaranno.csv) | is.null(opt$ctrlinf.filename)) {
  	cat("listvaranno_to_ctrlinf.R [OPTS] --listvaranno-csv <IN.csv> --ctrlinf-filename <OUT.ctrlinf>\n")
  	q()
  }

  anno <- read.csv(opt$listvaranno.csv, check.names=F, stringsAsFactors=F)

  # remove all variants on chrX, autosomal analysis only
  if (opt$include.chrX == F) {
    anno <- subset(anno, substr(anno[["Variant ID"]],1,2) != "X-")
  }

  # subset on annotations
  if (file.exists(opt$var.functions) ) { 
    opt$var.functions <- scan(opt$var.functions, what=character()) 
  } else {
    opt$var.functions <- strsplit(opt$var.functions, ",")[[1]]
  }
  anno <- subset(anno, Function %in% opt$var.functions)

  # subset on ppn humvar thresholds / classifs
  if (opt$min.ppn.humvar.score > 0) {
      anno <- anno[ is.na(anno[["Polyphen Humvar Score"]]) | 
                    anno[["Polyphen Humvar Score"]] > opt$min.ppn.humvar.score, ,drop=F]
  }
  if (is.null(opt$ppn.humvar.classifs) == F) {
      opt$ppn.humvar.classifs <- strsplit(opt$ppn.humvar.classifs, ",")[[1]]
      anno <- anno[ is.na(anno[["Polyphen Humvar Prediction"]]) | 
                    anno[["Polyphen Humvar Prediction"]] %in% opt$ppn.humvar.classifs, ,drop=F]
  }
  # subset on mean coverage in ExAC
  if (is.null(opt$exac.mean.cov.min) == F) {
	  anno <- anno[ anno[["ExAC Mean Coverage"]] >= opt$exac.mean.cov.min, , drop = F]
  }  
  
  # subset on MAF for every defined ExAC subgroup
  if (opt$maf.threshold < 1) {
    opt$maf.threshold.groups <- strsplit(opt$maf.threshold.groups, ",")[[1]]
    for (group in opt$maf.threshold.groups) {
      group.colname <- paste("ExAC", group, "maf", collapse="_")
  	  anno <- anno[ anno[[group.colname]] <= opt$maf.threshold, , drop = F]
    }
  }

  # get vector with counts of 0X/1X/2X allele
  maf.allele.freqs.group <- "global"
  maf.colname <- paste(c("ExAC", maf.allele.freqs.group, "maf"), collapse=" ")
  allele.count.colname <- paste(c("ExAC", maf.allele.freqs.group, "gts"), collapse=" ")  

  anno <- anno[ is.na( anno[[maf.colname]] ) == F, , drop = F]
  
  na.count.rows <- which( is.na(anno[[allele.count.colname]]) )
  anno[na.count.rows, allele.count.colname] <- paste(anno[na.count.rows, "ExAC Sample Covered 10x"], 
  													 "0", "0", sep="/")

  # get rid of unwanted single quotations
  anno[["Gene Name"]] <- gsub("'", "", anno[["Gene Name"]])
  anno[[allele.count.colname]] <- gsub("'", "", anno[[allele.count.colname]])  

  # build up ctrlinf data frame
  allele.counts.list <- strsplit(anno[[allele.count.colname]], "/")
  allele.0x.counts <- sapply(allele.counts.list, function(x) { return( as.numeric(x[1]) ) } )
  allele.1x.counts <- sapply(allele.counts.list, function(x) { return( as.numeric(x[2]) ) } )
  allele.2x.counts <- sapply(allele.counts.list, function(x) { return( as.numeric(x[3]) ) } )

  sample.recodeA.names <- paste(anno[["Variant ID"]], anno[["Alt Allele"]], sep="_")
  inv.af <- which(allele.0x.counts < allele.2x.counts)
  sample.recodeA.names[inv.af] <- paste(anno[inv.af, "Variant ID"], anno[inv.af, "Allele"], sep="_")
  
  if (opt$rvTDT.format == T) {
    controlinf <- data.frame(allele.2x.counts=allele.2x.counts,
                             allele.1x.counts=allele.1x.counts,
                             allele.0x.counts=allele.0x.counts,
                             ctrl.mean.cov=anno[["ExAC Mean Coverage"]],
                             gene=anno[["Gene Name"]])
    
  } else {
    controlinf <- data.frame(gene=anno[["Gene Name"]], varID=sample.recodeA.names,
  	    					   keep.var=TRUE, 
  		    				   allele.2x.counts=allele.2x.counts,
  			    			   allele.1x.counts=allele.1x.counts,  
  				        	   allele.0x.counts=allele.0x.counts,  
  						       ctrl.mean.cov=anno[["ExAC Mean Coverage"]])
  }					   
  rownames(controlinf) <- sample.recodeA.names
  
  write.table(controlinf, file=opt$ctrlinf.filename, row.names=T, col.names=F,
  			  quote=F, sep="\t")

}

GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--listvaranno-csv"), action="store", type="character",
                default=NULL, dest="listvaranno.csv", 
                help="name of required listvaranno.csv file, [default %default]"),

    make_option(c("--ctrlinf-filename"), action="store", type="character",
                default=NULL, dest="ctrlinf.filename", 
                help="name of required .ctrlinf output file, [default %default]"),

    make_option(c("--maf-threshold"), action="store", type="double",
                default=1, dest="maf.threshold", 
                help="MAF threshold in cohort data to use, [default %default]"),

    make_option(c("--var-functions"), action="store", type="character",
                default="START_LOST,STOP_GAINED,STOP_LOST,NON_SYNONYMOUS_CODING,NON_SYNONYMOUS_START,START_GAINED,SPLICE_SITE_ACCEPTOR,SPLICE_SITE_DONOR,FRAME_SHIFT,CODON_INSERTION,CODON_CHANGE_PLUS_CODON_INSERTION,CODON_DELETION,CODON_CHANGE_PLUS_CODON_DELETION",
                dest="var.functions", 
                help="annotated functions required for inclusion of a variant [default %default]"),

    make_option(c("--min-ppn-humvar-score"), action="store", type="double",
                default=0, dest="min.ppn.humvar.score",
                help="minimum polyphen humvar score required, [default %default]"),

    make_option(c("--ppn-humvar-classifs"), action="store", type="character",
                default=NULL, 
                dest="ppn.humvar.classifs",
                help="polyphen humvar classifications allowed, [default %default]"),

    make_option(c("--maf-threshold-groups"), action="store", type="character",
                default="global,afr,amr,eas,sas,fin,nfe", dest="maf.threshold.groups", 
                help="ExAC subgroups required to meet MAF threshold, [default %default]"),

    make_option(c("--maf-allele-freqs-group"), action="store", type="character",
                default="global", dest="maf.allele.freqs.group", 
                help="ExAC subgroup to use allele freqs for, [default %default]"),

   make_option(c("--exac-mean-cov-min"), action="store", type="integer",
                default=NULL, dest="exac.mean.cov.min", 
                help="minimum mean coverage for variant across ExAC samples,
                [default %default]"),

   make_option(c("--include-chrX"), action="store_true", default=FALSE,
               dest="include.chrX",
               help="include chrX, [default %default]"),

   make_option(c("--rvTDT-format"), action="store_true", default=FALSE,
               dest="rvTDT.format", 
               help="format output for rvTDT, [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
    main()
}
