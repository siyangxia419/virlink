library(virlink)

data("peptide_z")
peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)

test_that("output is a data frame", {
  test_result <- peptide_pairwise_cooccurrence(d = peptide_bi,
                                               test_method = "fisher",
                                               correlation = TRUE)

  expect_equal(nrow(test_result), ncol(peptide_bi) * (ncol(peptide_bi) - 1) / 2)
  expect_equal(unique(c(test_result$id1, test_result$id2)), as.numeric(colnames(peptide_bi)))
  expect_true(min(test_result$phi, na.rm = TRUE) >= -1 & max(test_result$phi, na.rm = TRUE) <= 1)
  expect_true(min(test_result$prop_1_in_2, na.rm = TRUE) >= -1 & max(test_result$prop_1_in_2, na.rm = TRUE) <= 1)
  expect_true(min(test_result$prop_2_in_1, na.rm = TRUE) >= -1 & max(test_result$prop_2_in_1, na.rm = TRUE) <= 1)
})
#> Test passed
