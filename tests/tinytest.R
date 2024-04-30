if (requireNamespace("tinytest", quietly=TRUE) ) {
  if (length(unclass(packageVersion("plmmr"))[[1]]) == 4 | Sys.getenv('R_FORCE_TEST') == 'TRUE') {
    tinytest::test_package("plmmr", pattern="^[^_].*\\.[rR]$")
  }
}
