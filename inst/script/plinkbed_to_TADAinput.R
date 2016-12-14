#!/usr/bin/Rscript

require(optparse)
require(data.table)
require(snpStats)

require(rvatk)

Main <- function() {
  option.list <- GetOptList()
  opt <- parse_args(OptionParser(option_list=option.list))

  req.args <- c("bed.file", "varinfo.file", "dnm.file",
                "gene.mu.file","trio.ped")
  for (arg.name in req.args) {
    if (is.null(opt[[arg.name]])) {
      cat("plinkbed_to_TADAinput.R 
              --bed-file <FILE> 
              --varinfo-file <FILE>
              --dnm-file <FILE>
              --gene-mu-file <FILE>
              --trio-ped <FILE>
              [--sib-ped FILE]
              [--casectrl-ped FILE]
              [--qc-proportion DOUBLE]\n")
      q()
    }
  }

  ### HANDLING OF VARINFO FILE FIRST
  ## read in varinfo file
  ## establish discrete var classes from varinfo file
  ## read in mutation rate file, establish that var classes sync with varinfo  

  ### READING IN OF LARGE FILES
  ## read in plink bed file, will end up reading bim+fam files by default

  ### HANDLING OF TRIO INFO 
  ## read in trio ped
  ## create a trio-only subset of plink bed file
  ## de novo mutation data handling
  # retrieve mutation rate per var class from mut rate table
  # retrieve mutation count per class from dnm file, varinfo table intersect
  # (store above info in list)
  ## inherited var data handling
  # for each gene,retrieve number of parental het vars for each class
  # count number of these vars transmitted to proband
  # subtract n transmitted from total parental vars to get ntrans count
  # store trans, ntrans vals for each gene/varclass combn, in a list
  
  ### HANDLING OF CASECTRL INFO
  ## read in casectrl ped
  ## create a casectrl-only subset of plink bed file
  ## for each gene/varclass combn, count total var in cases, in ctrls, store

  ### HANDLING OF SIB INFO
  ## read in sib ped, pheno==2 for proband, pheno==1 for affected sib
  ## create a sib-only subset of plink bed file
  ## for each gene/varclass combn, count total var in proband
  ## also count total var that are in proband + sib
  ## expected ratio of proband only vs. proband+sib ~ 1/1 
  
}


GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--bed-file"), action="store", type="character",
                default=NULL, dest="bed.file", 
                help="name of required plink bed file, [default %default]"),

    make_option(c("--varinfo-file"), action="store", type="character",
                default=NULL, dest="varinfo.file",
                help="name of required varinfo file, [default %default]"),

    make_option(c("--dnm-file"), action="store", type="character",      
                  default=NULL, dest="dnm.file",                          
                  help="name of required dnm file, [default %default]"), 

    make_option(c("--gene-mu-file"), action="store", type="character",
                  default=NULL, dest="gene.mu.file", 
                  help="name of required gene mu file, [default %default]"),

    make_option(c("--trio-ped"), action="store", type="character",
                default=NULL, dest="trio.ped",
                help="name of trio ped file, [default %default]"),

    make_option(c("--sib-ped"), action="store", type="character", 
                default=NULL, dest="sib.ped", 
                help="name of sibling ped file, [default %default]"), 

    make_option(c("--casectrl-ped"), action="store", type="character",
                default=NULL, dest="casectrl.ped",
                help="name of casectrl ped file, [default %default]"), 

    make_option(c("--qc-proportion"), action="store", type="double",
                default=0.8, dest="qc.proportion",
                help="proportion of cases required at each var 
                      for usage in TDT scenarios, [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
  Main()
}
