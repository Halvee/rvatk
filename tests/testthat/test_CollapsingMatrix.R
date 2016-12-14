
# initialize collapsing matrix to compare against
clpsmat <- data.frame("sample/gene"=c("GENE1","GENE2","GENE3"),
                      sample1=c(1,0,0),
                      sample5=c(1,1,0),
                      ctrl1=c(1,0,0),
                      ctrl2=c(0,1,0),
                      ctrl3=c(0,0,0),
                      ctrl4=c(0,1,0),
                      ctrl5=c(0,0,0),
                      ctrl6=c(0,1,0))
rownames(clpsmat) <- clpsmat[,1]
clpsmat[,1] <- NULL


context("CollapsingMatrix tests")

# test collapsing matrix read function
test_that("collapsing matrix read function works properly", {
  mat <- CollapsingMatrixRead("casectrl.clpsmat")
  expect_equal(clpsmat, mat)
})

# convert qualvar, use to initialize collapsing matrix 
test_that("collapsing matrix init. and read from qualvar works properly", {
  geneset <- c("GENE1","GENE2","GENE3")
  samples <- c("sample1","sample5","ctrl1","ctrl2","ctrl3","ctrl4","ctrl5","ctrl6")
  qualvar=data.frame("Gene Name"=c(rep("GENE1",5),rep("GENE2",6)),
                     "Variant ID"=c(rep("1-100-A-C",2), "1-105-C-T",
                                    rep("1-110-G-C",2), "2-100-A-C",
                                    "2-106-C-T", rep("2-111-G-A",3),
                                    "2-112-G-T"),
                     "Sample Name"=c("sample1","sample3","ctrl1","sample5",
                                     "sample8","ctrl6","ctrl4","sample5",
                                     "sample6","sample8","ctrl2"),
                     "Genotype"=c(rep(1,9),2,1))
  mat <- CollapsingMatrixInit(geneset, samples)
  mat <- CollapsingMatrixRecode(mat, qualvar)
  expect_equal(clpsmat, mat)
})
