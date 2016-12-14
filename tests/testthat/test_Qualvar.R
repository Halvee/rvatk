
context("Qualvar tests")

# test that file read capability of Qualvar
qualvar=data.frame("V1"=c(rep("GENE1",5),rep("GENE2",6)),
                   "V2"=c(rep("1-100-A-C",2), "1-105-C-T",
                          rep("1-110-G-C",2), "2-100-A-C",
                          "2-106-C-T", rep("2-111-G-A",3),
                          "2-112-G-T"),
                   "V3"=c("sample1","sample3","ctrl1","sample5",
                          "sample8","ctrl6","ctrl4","sample5",
                          "sample6","sample8","ctrl2"),
                   "V4"=c(rep(1,9),2,1),
                   stringsAsFactors=F, check.names=F)

test_that("Qualvar read function works properly", {
  qualvar.test <- QualvarRead("famctrl.qualvar")
  expect_equal(qualvar, qualvar.test)
})

# test sample / geneset subsetting capability of Qualvar
qualvar.casectrl<-data.frame("V1"=c(rep("GENE1",3),rep("GENE2",4)),
                             "V2"=c(rep("1-100-A-C",1), "1-105-C-T",
                                    rep("1-110-G-C",1), "2-100-A-C",
                                    "2-106-C-T", rep("2-111-G-A",1),
                                    "2-112-G-T"),
                             "V3"=c("sample1","ctrl1","sample5",
                                    "ctrl6","ctrl4","sample5",
                                    "ctrl2"),
                             "V4"=c(rep(1,7)),
                             stringsAsFactors=F, check.names=F)
test_that("Qualvar subsetting is working properly", {
  qualvar.test <- QualvarRead("famctrl.qualvar",
                               samples=c("sample1","sample5","ctrl1",
                                         "ctrl2","ctrl3",
                                         "ctrl4","ctrl5","ctrl6"),
                               geneset=c("GENE1","GENE2","GENE3") )
  rownames(qualvar.casectrl) <- NULL
  rownames(qualvar.test) <- NULL
  expect_equal(qualvar.casectrl, qualvar.test)
})
