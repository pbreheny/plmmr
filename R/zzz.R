.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is version ", packageVersion(pkgname),
                        " of ", pkgname, ".\n",
                        "Note: ", pkgname, " depends on the package bigalgebra, the current GitHub version of which is throwing some warnings for filebacked analysis.
\n https://github.com/fbertran/bigalgebra/issues/2 \n
If you see a warning about 'stack imbalance' appear while you are fitting a model with plmm() or cv_plmm(), we recommend downloading the last stable version of bigalebra.
If you are analyzing data stored in your R session memory, you can disregard this message.")
}