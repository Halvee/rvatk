
require(optparse)
require(rvatk)

main <- function() {
  
  # get opts, args
  option.list <- GetOptList()
  optarg <- parse_args(OptionParser(option_list=option.list),
                        positional_arguments = TRUE)
  opt <- optarg$options
  qv.files <- optarg$args

  # read in sample geneset
  geneset <- scan(opt$geneset.file, what=character())

  # read in covariate file, sample ID is expected to be column 1                
  covar <- read.table(opt$sample.covar, header=T, row.names=1)
  samplenames <- rownames(covar)

  # read each qualvar file, store df's in list, link with qualvar.names
  qv.sets <- list()
  qv.names <- strsplit(opt$qualvar.names,",")[[1]]
  stopifnot(length(qv.names) == length(qv.files))
  for (i in 1:length(qv.names)) {
    x <- QualvarRead(qv.files[i],samples=samplenames, geneset=geneset)
    qv.sets[[ qv.names[i] ]] <- QualvarRead(qv.files[i])
  }

  # read in covariate file, sample ID is expected to be column 1
  covar <- read.table(opt$sample.covar, header=T, row.names=1)

  # take casectrl column and make sure scaling is 0/1, not 1/2
  for (binary.pred in c("CASECTRL","SEX")) {
    stopifnot( length(intersect(c(binary.pred), 
                                colnames(covar))) == 1)
    if (max(covar[[binary.pred]]) == 2) {
      covar[[binary.pred]] <- covar[[binary.pred]] - 1
    }
  }

  # add columns of qual var counts per sample
  covar$qv.all <- 0
  for (qv.name in qv.names) {
    covar[[paste("qv",qv.name,sep=".")]] <- 0
  }

  # get sum of qualifying variants per sample across all qv files, and also
  # create a qv df of all unique var calls across qv files
  qv.all <- c()
  for (qv.name in qv.names) {
    print(qv.name)
    qv <- qv.sets[[ qv.name ]]
    qv.all <- rbind(qv.all, qv.sets[[qv.name]])
    qv.all <- unique(qv.all)
  }

  # store the total number of unique var calls per sample across qv files
  for (sample.name in unique(qv.all[, 3])) {
    covar[sample.name, "qv.all"] <- sum(qv.all[ qv.all[,3] == sample.name, 4])
  }

  # get number of vars per qv set, subsetted on geneset if user provided one
  geneset <- c()
  if (is.null(opt$geneset.file) == F) {
    geneset <- scan(opt$geneset.file, what=character())
  }
  for (qv.name in qv.names) {
    qv <- qv.sets[[ qv.name ]]
    if (length(geneset) > 0) {
      qv <- qv[ qv[,1] %in% geneset, ]
    }
    covar.col <- paste("qv",qv.name,sep=".")
    for (sample.name in unique(qv[, 3])) {
      covar[sample.name, covar.col] <- sum(qv[ qv[,3] == sample.name, 4])
    }
  }

  # if maximum variant count defined, filter samples from analysis that exceed threshold
  if (is.null(opt$max.sample.varcount) == F) {
    print(covar[covar$qv.all > opt$max.sample.varcount, ])
    covar <- covar[ covar$qv.all <= opt$max.sample.varcount, ]
  }
  
  # load expression for statistical model we're using
  stopifnot(is.null(opt$mdl.eqn) == F)
  if (file.exists(opt$mdl.eqn) == T) {
    expr.str <- scan(opt$mdl.eqn, what=character(), sep="\n")[1]
  } else {
    expr.str <- opt$mdl.eqn
  }

  # initialize base logistic regression model
  expr.str <- paste(expr.str, "qv.all", sep=" + ")
  print(expr.str)
  expr <- as.formula(expr.str)
  mdl.null <- glm(expr, data = covar, family="binomial")
  for (qv.name in qv.names) {
    expr.str.alt <- paste(expr.str,paste("qv",qv.name,sep="."), sep=" + ")
    print(expr.str.alt)
    expr.alt <- as.formula(expr.str.alt)
    mdl.alt <- glm(expr.alt, data = covar, family="binomial")
    #res <- lm(expr.alt, data = covar)
    odds.ratios <- exp(coef(mdl.alt))
    print(odds.ratios)
    print(summary(mdl.alt))
    x<-anova(mdl.null, mdl.alt)
    print(summary(x))
  }

}

GetOptList <- function() {
  # Get option list.
  #
  # Returns:
  #   Option list to parse options.
  option.list <- list(
    make_option(c("--qualvar-names"), action="store", type="character",
                default=NULL, dest="qualvar.names", 
                help="names assigned to each qualvar set, comma delimited, 
                      [default %default], required"),

    make_option(c("--sample-covar"), action="store", type="character",
                default=NULL, dest="sample.covar",
                help="file with covariates per sample,
                      [default %default], required"),

    make_option(c("--geneset-file"), action="store", type="character", 
                default=NULL, dest="geneset.file",
                help="name of geneset file, [ddefault %default], required"),

    make_option(c("--max-sample-varcount"), action="store", type="integer",
                default=NULL, dest="max.sample.varcount",
                help="maximum number of variants allowed for a sample to
                      be included in analysis, [default %default]"),

    make_option(c("--mdl-eqn"), action="store", type="character",
                default=NULL, dest="mdl.eqn",
                help="equation to use for stat model, required, [default %default]")

  )

  return(option.list)
}

if (interactive() == F) {
  main()
}
