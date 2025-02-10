.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is ", pkgname, " ", utils::packageVersion(pkgname), ".\n")
}

# Define an environment to store the state
.plmmr_env <- new.env()

# Initialize the state
.plmmr_env$warning_shown <- FALSE