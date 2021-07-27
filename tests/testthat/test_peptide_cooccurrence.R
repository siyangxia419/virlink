library(virlink)

data("peptide_z")
peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)

test_that("output is a numeric vector", {
  test_result1 <- peptide_cooccurrence(x1 = peptide_bi[, 1],
                                       x2 = peptide_bi[, 2],
                                       test_method = "fisher",
                                       correlation = TRUE)

  expect_equal(length(test_result1), 7)
  expect_equal(names(test_result1), c("TT", "TF", "FT", "FF", "oddsr", "p", "phi"))
  expect_equal(sum(test_result1[c("TT", "TF", "FT", "FF")]), length(peptide_bi[, 1]))

  test_result2 <- peptide_cooccurrence(x1 = peptide_bi[, 1],
                                       x2 = peptide_bi[, 2],
                                       test_method = "chisq",
                                       correlation = FALSE)

  expect_equal(length(test_result2), 6)
  expect_equal(names(test_result2), c("TT", "TF", "FT", "FF", "xsq", "p"))
  expect_equal(sum(test_result2[c("TT", "TF", "FT", "FF")]), length(peptide_bi[, 2]))
})
#> Test passed
