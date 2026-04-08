# Oscar Rysavy - tests for various data functions

# Test 1 - Example datasets can be unzipped -----------------------
tmp <- tempdir()
expect_message(unzip_example_data(tmp), pattern = "Unzipped")
