# GLOBAL VARS
kGeneNameCol <- "Gene Name"
kVariantIdCol<- "Variant ID"
kSampleNameCol <- "Sample Name"
kGenotypeCol <- "Genotype"
kLvgQualvarCols <- c(kGeneNameCol, kVariantIdCol, kSampleNameCol)
kCphtGeneNameCol <- "Gene"
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


#' read in a file where each line is a seperate entry into a vector.
#'
#' @param filename name of input file.
#' 
#' @export
ReadListFile <- function(filename) {
  mylist <- scan(filename, what=character(), quiet=T)
  return(mylist)
}

#' read in a .condition file, where each line represents conditions to apply.
#'
#' @param condition.filename name of .condition file.
#' @export
ReadConditionFile <- function(condition.filename) {
  cnd <- read.table(condition.filename, stringsAsFactors=F)
  colnames(cnd) <- c("col.name", "rel.operator", "val")
  return(cnd)
}

#' add a leave-one-out minor allele freq. column with values to genotype table.
#'
#' @param geno data frame of genotypes.
#' @param n.samples total number of samples.
#' @export
AddLooMafCol <- function(geno, n.samples, variant.id.col="Variant ID") {
  geno[["loo maf"]] <- 0
  var.loo.counts <- table(geno[[variant.id.col]]) - 1
  var.loo.counts <- var.loo.counts[var.loo.counts > 0]
  if (length(var.loo.counts) > 0) {
    unique.counts <- sort(as.numeric(unique(var.loo.counts)))
    max.loo.count <- max(var.loo.counts)
    for (i in unique.counts) {
      vars <- names(var.loo.counts[var.loo.counts == i])
      loo.maf.i <- as.numeric(i / n.samples)
      geno[ geno[[variant.id.col]] %in% vars,
            "loo maf"] <- loo.maf.i
    }
  }
  return(geno)
}

#' subset on table based on rows in condition data frame.
#'
#' @param tbl data frame of values with column names.
#' @param cnd data frame of conditions to subset tbl on.
#' @export
TblCndSubset <- function(tbl, cnd) {
  for (i in 1:nrow(cnd)) {
    col.name <- cnd[i,"col.name"]
    rel.operator <- cnd[i,"rel.operator"]
    val <- cnd[i,"val"]
    tbl <- ApplyThreshold(tbl, col.name, rel.operator, val)
  }
  return(tbl)
}

#' split string to vector
#' 
#' @param str input string.
#' @param type type of variable that each output vector element should be.
#' @param delim delimiter character in input string.
#' @export
StringToSet <- function(str,
                        type="character",
                        delim=",") {
  if (is.null(str)==F) {
    str.set <- strsplit(str, delim)[[1]]
  } else {
    str.set <- c()
  }
  if (type!="character") {
    for (i in 1:length(str.set)) {
      str.set[i] <- as(str.set[i], type)
    }
  }
  return(str.set)
}

#' subset on input table on a user defined condition.
#' 
#' @param tbl input data frame.
#' @param col.name column name to subset input data frame on.
#' @param cmp string representing which operator to apply on df column.
#' @param val value linked to input tbl by cmp operator.
#' @export
ApplyThreshold <- function(tbl, col.name, cmp, val,
                           na.numeric.convert=0) {
  if (cmp == "in") {
    valset <- StringToSet(val)
    tbl <- tbl[ which(tbl[[col.name]] %in% valset), , drop=F]
  } else if (cmp == "notin") {
    valset <- StringToSet(val)
    tbl <- tbl[ which((tbl[[col.name]] %in% valset)==F), , drop=F]
  } else if (cmp == "grep") {
    tbl <- tbl[ which(grepl(val, tbl[[col.name]]) == T), , drop=F]
  } else if (cmp == "grepv") {
    tbl <- tbl[ which(grepl(val, tbl[[col.name]]) == F), , drop=F]
  } else if (cmp == "eq" | cmp == "==") {
    tbl <- tbl[ which(tbl[[col.name]] == val), , drop =F]
  } else if (cmp == "noteq" | cmp == "!=") {
    tbl <- tbl[ which(tbl[[col.name]] != val), , drop =F]
  } else if (cmp == "gt" | cmp == ">") {
    val <- as.double(val)
    valset <- as.double( tbl[[col.name]] )
    valset <- ifelse(is.na(valset), na.numeric.convert, valset)
    tbl <- tbl[ which(valset > val), , drop =F]
  } else if (cmp == "gte" | cmp == ">=") {
    val <- as.double(val)
    valset <- as.double( tbl[[col.name]] )
    valset <- ifelse(is.na(valset), na.numeric.convert, valset) 
    tbl <- tbl[ which(valset >= val), , drop =F]
  } else if (cmp == "lt" | cmp == "<") {
    val <- as.double(val)
    valset <- as.double( tbl[[col.name]] )
    valset <- ifelse(is.na(valset), na.numeric.convert, valset) 
    tbl <- tbl[ which(valset < val), , drop =F]
  } else if (cmp == "lte" | cmp == "<=") {
    val <- as.double(val)
    valset <- as.double( tbl[[col.name]] )
    valset <- ifelse(is.na(valset), na.numeric.convert, valset) 
    tbl <- tbl[ which(valset <= val), , drop =F]
  }
  return(tbl)
}

#' read a large table (can be gzipped and/or GB in size) into a data frame.
#' @param filename name of table file.
#' @param ... args to be passed to data.table fread function
#' @importFrom data.table fread
#' @export
ReadLargeTable <- function(filename, ...) {
  if (grepl(".gz$",filename)==T) {
    filename.full <- paste("gunzip -c",
                      filename,
                      sep=" ")
  } else {
    filename.full <- filename
  }
  mat <- fread(filename.full,
               data.table=F,
               ...)
  return(mat)
}

#' write large table to output file, with support for gzip.
#' 
#' @param tbl data frame to be written to file.
#' @param filename name of output file to write to, ".gz" at end for gzip.
#' @param ... parameters to be passed on to write.table function.
#' @export
WriteLargeTable <- function(tbl, filename, ...) {
  if (grepl(".gz$", filename) == T) {
    filename.full <- gzfile(filename)
  } else {
    filename.full <- filename
  }
  write.table(tbl, file=filename.full,
              ...)
}

#' perform a case/control burden test on aggregated rare variant sequence data.
#' 
#' @param mat collapsing matrix data structure.
#' @param genes list of genes to include in burden test.
#' @param cases vector of case names.
#' @param ctrls vector of ctrl names.
#' @param collapse.thresh integer threshold for inclusion of sample genotype,
#'                        if >0 then collapsing analyses are used, else 
#'                        poisson tests used (case vs. ctrl variant rate).
#' @param ... additional args for fisher.test or poisson.test
#' @export
BurdenTest <- function(mat, genes, 
                       cases, ctrls,
                       collapse.thresh=0, 
                       ...) {
  n.cases <- length(cases)
  n.ctrls <- length(ctrls)
  cases <- intersect(cases, colnames(mat))
  ctrls <- intersect(ctrls, colnames(mat))
  genes <- intersect(genes, rownames(mat))

  res <- list()
  mat.g.counts <- colSums(mat[genes, , drop=F])
  if (collapse.thresh == 0) {
    mat.g.counts.cases <- sum(mat.g.counts[cases])
    mat.g.counts.ctrls <- sum(mat.g.counts[ctrls])
    var.rate.cases <- mat.g.counts.cases / n.cases
    var.rate.ctrls <- mat.g.counts.ctrls / n.ctrls
    res.full <- poisson.test(mat.g.counts.cases, n.cases,
                             var.rate.ctrls, ...)
    res$n.var.case <- mat.g.counts.cases
    res$n.case <- n.cases
    res$var.rate.case <- var.rate.cases
    res$n.var.ctrl <- mat.g.counts.ctrls
    res$n.ctrl <- n.ctrls
    res$var.rate.ctrl <- var.rate.ctrls
    res$conf.int.lower <- res.full$conf.int[1]
    res$conf.int.upper <- res.full$conf.int[2]
    res$p.value <- res.full$p.value

  } else {
    mat.g.counts <- ifelse(mat.g.counts >= collapse.thresh, 1, 0)
    mat.g.counts.cases <- sum(mat.g.counts[cases])
    mat.g.counts.ctrls <- sum(mat.g.counts[ctrls])
    n.case.qual <- mat.g.counts.cases
    n.case.nonqual <- n.cases - mat.g.counts.cases
    n.ctrl.qual <- mat.g.counts.ctrls
    n.ctrl.nonqual <- n.ctrls - mat.g.counts.ctrls
    FET.tbl <- data.frame(case=c(n.case.qual, n.case.nonqual),
                          ctrl=c(n.ctrl.qual, n.ctrl.nonqual))
    res.full <- fisher.test(FET.tbl, ...)
    res$n.case.qual <- n.case.qual
    res$n.case.nonqual <- n.case.nonqual
    res$n.ctrl.qual <- n.ctrl.qual
    res$n.ctrl.nonqual <- n.ctrl.nonqual
    res$odds.ratio <- res.full$estimate
    res$conf.int.lower <- res.full$conf.int[1]
    res$conf.int.upper <- res.full$conf.int[2]
    res$p.value <- res.full$p.value
  }
  res$method <- ProcessString(res.full$method)
  res$alternative <- ProcessString(res.full$alternative)

  return(res)
}

#' convert process input string into a form more agreeable with R.
#' @param str input string.
#' @export
ProcessString <- function(str) {
  str <- gsub(" ",".",str)
  str <- gsub("'","",str)
  return(str)
}

#' read genesets file to a list of genesets.
#' @param genesets.file name of .genesets file.
#' @export
ReadGenesetsFile <- function(genesets.file) {
  genesets <- list()
  fh <- file(genesets.file, open="r")
  lines <- readLines(fh)
  rows <- c()
  for (i in 1:length(lines)){
    row.i <- strsplit( lines[i], "\t" )[[1]]
    geneset.name <- row.i[1]
    geneset <- strsplit(row.i[2], "," )[[1]]
    genesets[[geneset.name]] <- geneset
  }
  close(fh)
  return(genesets)
}

#' update input table with values of res.list at row i
#' @param tbl input table.
#' @param row.i row number in table to edit.
#' @param res.list list of results to iterate through and add vals to table.
#' @export
UpdateTable <- function(tbl, row.i, res.list) {
  for (col.name.j in colnames(tbl)) {
    tbl[row.i, col.name.j] <- res.list[[col.name.j]]
  }
  return(tbl)
}

#' table sampleped data frame and get list of individual ID -> family ID.
#' @param sampleped data frame of sample data in sampleped format.
#' @export
SamplepedIidFidList <- function(sampleped) {
  iid.fid <- list()
  for (i in 1:nrow(sampleped)) {
    iid.fid[[ sampleped[i, 2] ]] <- sampleped[i,1]
  }
  return(iid.fid)
}

#' for input variant calls and genotypes, return genotype counts, taking into
#' account the chromosome and sample gender for each variant call.   For 
#' example, a 'homozygous' on an x chromosome of a male would be an allele 
#' count of only 1, rather than 2, whereas a female homozygous call on an x
#' chromosome would be an allele count of 2.Such rules do not apply to 
#' autosomes, where each sample should have two copies of each.
#' @param sample.ID vector of samples IDs
#' @param variant.ID vector of variant IDs
#' @param genotype vector of genotype strings (het / hom)
#' @param sampleped data frame of samples in sampleped format
#' @export
GetAlleleCounts <- function(sample.ID, 
                            variant.ID, 
                            genotype,
                            sampleped) {
  genotype <- ifelse(genotype=="hom",2,1)
  calls <- data.frame(sample.ID=sample.ID,
                      variant.ID=variant.ID,
                      genotype=genotype,
                      stringsAsFactors=F)
  samplegender <- list()
  for (i in 1:nrow(sampleped)) {
    sample.ID.i <- sampleped[i,2]
    sample.gender.i <- as.numeric(sampleped[i,5])
    samplegender[[ sample.ID.i ]] <- sample.gender.i
  }
  gender <- as.numeric(unlist(samplegender[ calls$sample.ID ]))
  calls$gender <- gender
  if ( nrow(calls[ (grepl("X-", calls$variant.ID)) &
                   (calls$gender == 1) & 
                   (calls$genotype == 2), ]) > 0) {
    calls[ grepl("X-", calls$variant.ID) &
           calls$gender == 1 &
           calls$genotype == 2, "genotype"] <- 1
  }
  return(calls$genotype)
}

LoadExpression<-function(expression.str, is.file=F) {
  if (is.file == T) {
    expr <- scan(expression.str, what=character(),sep="\n")
    expr <- expr[1]
  } else {
    expr <- expression.str
  }
  expr <- as.expression(expr)
  return(expr)
}
