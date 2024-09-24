#' A function to rotate filebacked data
#'
#' @returns a list with 4 items:
#'  * stdrot_X: `X` on the rotated and re-standardized scale
#'  * rot_y: `y` on the rotated scale (a numeric vector)
#'  * stdrot_X_center: numeric vector of values used to center `rot_X`
#'  * stdrot_X_scale: numeric vector of values used to scale `rot_X`
#'
#' @keywords internal
#'
rotate_filebacked <- function(prep, ...){
  w <- (prep$eta * prep$s + (1 - prep$eta))^(-1/2)
  wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")

  # rotate X
  std_X <- prep$std_X
  rot_X <- wUt%*%std_X # using %*% method from bigalgebra
  stdrot_X <- bigmemory::big.matrix(nrow = nrow(wUt), ncol = ncol(std_X))

  # rotate y
  rot_y <- wUt%*%prep$centered_y

  # re-standardize (since std_X is big, we do this in C++)
  std_rot <- .Call("big_std",
                   rot_X@address,
                   as.integer(count_cores()),
                   PACKAGE = "plmmr")
  stdrot_X@address <- std_rot[[1]]

  return(list(stdrot_X = stdrot_X,
              rot_y = rot_y,
              stdrot_X_center = std_rot[[2]],
              stdrot_X_scale = std_rot[[3]]))
}
