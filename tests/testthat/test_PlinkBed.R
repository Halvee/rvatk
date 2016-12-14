
context("PlinkBed tests")

# initialize PlinkBed obj
pb <- PlinkBedInit("famctrl.bed")

# get individual genotype calls, test by comparing against qualvar file
samples =c("sample1","sample3","ctrl1","sample5", 
           "sample8","ctrl6","ctrl4","sample5",
           "sample6","sample8","ctrl2")
vars = c(rep("1-100-A-C",2), "1-105-C-T",
         rep("1-110-G-C",2), "2-100-A-C",
         "2-106-C-T", rep("2-111-G-A",3),
         "2-112-G-T")
qualvar=data.frame("V1"=c(rep("GENE1",5),rep("GENE2",6)),
                   "V2"=vars,
                   "V3"=samples,
                   "V4"=c(rep(1,9),2,1),
                   stringsAsFactors=F, check.names=F)  
test_that("Qualvar genotype call retrieval is working properly", {
  samples.sub <- c("sample1","ctrl1","sample5","sample8")
  vars.sub <- c("1-100-A-C", "1-105-C-T", "2-111-G-A")
  pb.calls.gt1 <- PlinkBedCalls(pb, samples.sub, vars.sub, gt=1, as.df=T)
  pb.calls.gt2 <- PlinkBedCalls(pb, samples.sub, vars.sub, gt=2, as.df=T)
  pb.calls.gt1.exp <- data.frame(
                                 "X1"=c("1-100-A-C","1-105-C-T","2-111-G-A"),
                                 "X2"=c("sample1","ctrl1","sample5"),
                                 stringsAsFactors=F)
  pb.calls.gt2.exp <- data.frame("X1"=c("2-111-G-A"),
                                 "X2"=c("sample8"), 
                                 stringsAsFactors=F)
  colnames(pb.calls.gt1) <- NULL
  colnames(pb.calls.gt1.exp) <- NULL
  colnames(pb.calls.gt2) <- NULL
  colnames(pb.calls.gt2.exp) <- NULL
  expect_equal(pb.calls.gt1, pb.calls.gt1.exp)
  expect_equal(pb.calls.gt2, pb.calls.gt2.exp)
})
