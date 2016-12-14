

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
