library(biglasso)
data(colon)

# create toy data, add in an ID column and a simulated 'sex' variable
X <- colon$X
sex <- sample(1:2, nrow(X), replace = TRUE)
colon2 <- cbind(sex, X)
colnames(colon2)[1] <- c("sex")

# table of input data
write.table(colon2,
            "inst/extdata/colon2.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

# table with outcome
write.table(data.frame(ID = 1:length(colon$y),
                       y = colon$y),
            "inst/extdata/colon2_outcome.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)
