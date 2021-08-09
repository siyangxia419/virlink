library(virlink)

data("peptide_z")
peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)

test_that("output is a number if only phi coefficient is desired", {
  test_result1 <- peptide_cooccurrence(x = peptide_bi[, 1],
                                       y = peptide_bi[, 2],
                                       coefficient = "phi",
                                       perform_test = FALSE)
  expect_type(test_result1, type = "double")
})
#> Test passed

test_that("output is a numeric vector, otherwise", {
  test_result2 <- peptide_cooccurrence(x = peptide_bi[, 1],
                                       y = peptide_bi[, 2],
                                       coefficient = "both",
                                       perform_test = TRUE,
                                       test_method = "fisher")

  expect_type(test_result2, type = "double")
  expect_vector(test_result2)
  expect_equal(length(test_result2), 8)
  expect_equal(names(test_result2),
               c("TT", "TF", "FT", "FF", "oddsr", "p", "phi", "mean_prop"))
  expect_equal(sum(test_result2[c("TT", "TF", "FT", "FF")]),
               length(peptide_bi[, 1]))
})
#> Test passed

test_that("invalid input causes error", {
  expect_error(peptide_cooccurrence(x = peptide_bi[, 1], y = numeric(0)))
})
#> Test passed
