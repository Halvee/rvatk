

context("misc functions")

test_that("string processing", {
  teststrtoset <- "GENE1,GENE2,GENE3"
  expect_equal(StringToSet(teststrtoset), 
               c("GENE1","GENE2","GENE3"))
  teststrprocess <- c("'GENE 1'","'GENE 2'","'GENE 3'")
  expect_equal(ProcessString(teststrprocess),
               c("GENE.1","GENE.2","GENE.3"))
})

clpsmat <- data.frame("sample/gene"=c("GENE1","GENE2","GENE3"),
                      sample1=c(1,0,0),
                      sample5=c(1,1,0),
                      ctrl1=c(1,0,0),
                      ctrl2=c(0,1,0),
                      ctrl3=c(0,0,0),
                      ctrl4=c(0,1,0),
                      ctrl5=c(0,0,0),
                      ctrl6=c(0,1,0),
                      check.names=F,stringsAsFactors=F)

test_that("read large files", {
  clpsmat.test <- ReadLargeTable("casectrl.clpsmat",
                                 header=T, sep="\t",
                                 check.names=F,
                                 stringsAsFactors=F)
  expect_equal(clpsmat, clpsmat.test)
})
test_that("read large gzipped files", {
  clpsmat.test <- ReadLargeTable("casectrl.clpsmat.gz",
                                 header=T, sep="\t",
                                 check.names=F,
                                 stringsAsFactors=F)
  expect_equal(clpsmat, clpsmat.test)
})

test_that("read geneset files", {
  geneset <- c("GENE1","GENE2","GENE3")
  geneset.test <- ReadListFile("genes.geneset")
  expect_equal(geneset, geneset.test)
})

test_that("read .genesets files", {
  genesets <- list("genesetA"=c("GENE1","GENE2","GENE3"),
                   "genesetB"=c("GENE1","GENE4","GENE5"))
  genesets.test <- ReadGenesetsFile("genes.genesets")
  expect_equal(genesets, genesets.test)
})

test_that("condition file reading", {
  cnd <- data.frame(col.name=c("sampleID", 
                               "sampleID",
                               "varA", 
                               "varA"),
                    rel.operator=c("eq", "==",
                                   "lte", "<="),
                    val=c("sample1", "sample1",
                           "2", "2"), stringsAsFactors=F)
  cnd.file <- ReadConditionFile("test.condition")
  expect_equal(cnd, cnd.file)
})

test_that("subsetting based on condition tables (string subsetting)", {
  data <- data.frame(sampleID=paste("sample", seq(1,5), sep=""),
                     varA=c(1,2,3,1,2),
                     varB=c("X","Y","X","X","Y"),
                     varC=c(10.0, 10.1, 10.2,
                            10.1, 9.9), stringsAsFactors=F )
  rownames(data) <- NULL

  cnd.1 <- data.frame(col.name=c("varB"), 
                      rel.operator=c("=="), 
                      val=c("Y"), stringsAsFactors=F)
  symbols <- c("==", "!=")
  strs <- c("eq", "noteq")
  nrows <- c(2, 3)
  for (i in 1:2) {
    cnd.1[1, "rel.operator"] <- symbols[i]
    data.1.test <- TblCndSubset(data, cnd.1)
    expect_equal(nrow(data.1.test), nrows[i])
    cnd.1[1, "rel.operator"] <- strs[i]
    data.1.test <- TblCndSubset(data, cnd.1)
    expect_equal(nrow(data.1.test), nrows[i])
  }
  
  strs <- c("in", "notin")
  nrows <- c(2, 3)
  for (i in 1:2) {
    cnd.1[1, "rel.operator"] <- strs[i]
    data.1.test <- TblCndSubset(data, cnd.1)
    expect_equal(nrow(data.1.test), nrows[i])
  }

  strs <- c("grep","grepv")
  nrows <- c(2, 3)
  for (i in 1:2) {
    cnd.1[1, "rel.operator"] <- strs[i] 
    data.1.test <- TblCndSubset(data, cnd.1)
    expect_equal(nrow(data.1.test), nrows[i])
  }
})

test_that("subsetting based on condition tables (numeric subsetting)", {
  data <- data.frame(sampleID=paste("sample", seq(1,5), sep=""),
                     varA=c(1,2,3,1,2),
                     varB=c("X","Y","X","X","Y"),
                     varC=c(10.0, 10.1, 10.2,
                            10.1, 9.9), stringsAsFactors=F )
  rownames(data) <- NULL

  cnd.2 <- data.frame(col.name=c("varC"),
                      rel.operator=c("<"), 
                      val=c("10.1"), stringsAsFactors=F)
  symbols <- c("<", "<=", "==", "!=", ">=", ">")
  strs <- c("lt", "lte", "eq", "noteq", "gte", "gt")
  nrows <- c(2, 4, 2, 3, 3, 1)
  for (i in 1:5) {
    cnd.2[1, "rel.operator"] <- symbols[i]
    data.2.test <- TblCndSubset(data, cnd.2)
    expect_equal(nrow(data.2.test), nrows[i] )
    cnd.2[1, "rel.operator"] <- strs[i]
    data.2.test <- TblCndSubset(data, cnd.2)
    expect_equal(nrow(data.2.test), nrows[i] )
  }
})
