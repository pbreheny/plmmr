
### function to calculate population variance
### used in conjunction with our standardized X matrix
### since we assume that is our reference population and not a sample
varp <- function(x) mean((x-mean(x))^2)