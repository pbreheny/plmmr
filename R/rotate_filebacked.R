#' A function to rotate filebacked data 
#'
#'@returns a list with 4 items: 
#'  * rot_X" `X` on the rotated scale, without re-standardizing (an FBM)
#'  * rot_y: `y` on the rotated scale (a numeric vector)
#'  * stdrot_X: `X` on the rotated scale, re-standardized (an FBM)
#'  * stdrot_X_scale: numeric vector of values used to standardize `rot_X`
#'  
#' @keywords internal
#' 
rotate_filebacked <- function(prep, ...){
 browser()
  w <- (prep$eta * prep$s + (1 - prep$eta))^(-1/2)
  wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")

  # rotate X (since std_X is big, we do this in C++)
  std_X <- fbm2bm(prep$std_X)
  rot_X <- bigstatsr::FBM(nrow = nrow(wUt), ncol = ncol(std_X)) |> fbm2bm()
  rot_X_res <- .Call("rotate_filebacked",
                 std_X@address,
                 wUt, 
                 rot_X@address, 
                 as.integer(bigstatsr::nb_cores()),
                 PACKAGE = 'plmmr')
  rot_X@address <- rot_X_res[[1]]
  
  # rotate y
  rot_y <- wUt%*%prep$y
  
  # Pick up here: decide how to handle intercept and re-scaling
  
  # re-scale rot_X
  rot_X_scale_info <- bigstatsr::big_scale()(rot_X)
  
  # stdrot_X_center <- rot_X_scale_info$center
  stdrot_X_scale <- rot_X_scale_info$scale[-1] # NB: here, take off 'scale' corresponding to intercept column
  stdrot_X <- big_std(X = rot_X,
                      # center = stdrot_X_center, # NB: do not re-center rotated data
                      scale = stdrot_X_scale,
                      ns = 1:ncol(rot_X)-1) # NB: the -1 -> we took off the intercept
  
  return(list(rot_X = rot_X,
              rot_y = rot_y,
              stdrot_X = stdrot_X,
              stdrot_X_scale = stdrot_X_scale))
}