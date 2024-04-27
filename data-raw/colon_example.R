library(biglasso)
data(colon)

# create mock data with missing values 
sample(1:ncol(colon$X), 20) -> to_make_na
X <- colon$X
X[sample(1:nrow(X), 20), to_make_na] <- NA_real_
sum(is.na(X))

# create toy data, add in a simulated 'sex' variable 
sex <- sample(1:2, nrow(X), replace = T)
colon2 <- cbind(sex, colon$y, X)

write.table(colon2, 
            "inst/extdata/colon2.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
