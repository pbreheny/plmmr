# TKP
# August 2024
# Goal: test the big_cbind() function to assess what is causing data to be read
#   into memory. This function has been memory-intensive in test runs of plmm(),
#   even for filebacked data

n <- 1000
p <- 100000

# simulate data
A <- matrix(rnorm(n*3), n, 3) # A could be the matrix of unpenalized predictors from an external file (e.g., sex and age)
B <-  matrix(rnorm(n*p), n, p) |> bigmemory::as.big.matrix(backingfile = "test.bk",
                                                           type = 'double',
                                             backingpath = "inst/extdata/",
                                             descriptorfile = "test.desc")# B could be a matrix of genotype data

C <- bigmemory::filebacked.big.matrix(nrow = n,
                                      ncol = ncol(A) + p,
                                      type = 'double',
                                      backingfile = "combined_test.bk",
                                      backingpath = "inst/extdata/",
                                      descriptorfile = "combined_test.desc")

Rprof(tf <- "rprof.log", memory.profiling=TRUE)

# combine data
D <- big_cbind(A, B, C, quiet = F)

Rprof(NULL)

summaryRprof(tf)
# notice, the total RAM in the R session has increased... why is this?
# this can become quite challenging for GWAS data


# cleanup
file.remove(file.path("inst", "extdata",c("test.bk", "test.desc", "combined_test.bk", "combined_test.desc")))
