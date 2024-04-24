#' A function to convert the output of `process_plink()` into a `big.matrix` (such as could be passed to `biglasso::biglasso()`)
#'
#' @param path The file path to the RDS object containing the processed data. Do not add the '.rds' extension to the path. 
#' @param standardize Logical: should the data be returned in a column-standardized form? Defaults to TRUE. Please don't change this unless you are super confident about what you are doing...
#' @return A list with elements: 
#' * 'X': a `big.matrix` object with the design matrix of the data. If standardize = TRUE, the data corresponding to `std_X` is returned. 
#' Otherwise, the data corresponding to `subset_X` is returned 
#' * 'map': data frame with `.bim` file information 
#' * 'fam': data frame with `.fam` file information 
#' *  'ns': A vector indicating the which columns of X contain nonsingular features (i.e., features with variance != 0. 
#' * 'center': If standardize = TRUE, a vector of values for centering each column in X is returned 
#' * 'scale': If standardize = TRUE, a vector of values for scaling each column in X is returned 
#' @export
#'
#' @details This function is a wrapper combining `get_data()` and `fbm2bm()`
#'
#' @examples
#' \dontrun{
#' process_plink(data_dir = plink_example(parent = T),
#'   prefix = "penncath_lite",
#'   gz = TRUE,
#'   outfile = "process_penncath",
#'   # overwrite = TRUE, # uncomment if needed 
#'   impute_method = "mode")
#'   
#'   my_path <- paste0(plink_example(parent = T), "/penncath_lite")
#'   bm_data <- get_bm_from_plink(my_path)
#' }
get_bm_from_plink <- function(path, standardize = TRUE){
  
  rds <- get_data(path = path, returnX = FALSE)
  
  ret <- list(
    fam = rds$fam,
    map = rds$map,
    ns = rds$ns)
     
  if (standardize){
    ret$X <- fbm2bm(rds$std_X)
    ret$center <- rds$std_X_center
    ret$scale <- rds$std_X_scale
  } else {
    ret$X <- fbm2bm(rds$subset_X)
  }
  
  return(ret)
}