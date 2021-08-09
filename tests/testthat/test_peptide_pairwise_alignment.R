library(virlink)

data(peptide_df)

test_that("input is a data frame with the first two columns being 'id' and 'pep_aa'", {
  expect_error(peptide_pairwise_alignment(peptides = peptide_df[, -1]))
  expect_error(peptide_pairwise_alignment(peptides = peptide_df[, -2]))
  expect_error(peptide_pairwise_alignment(peptides = "NNNNNNN"))
})
#> Test passed
test_that("output has the right dimensions and structure", {
  test_result1 <- peptide_pairwise_alignment(peptides = peptide_df,
                                             self_comparison = TRUE,
                                             output_str = "data.table")

  expect_equal(nrow(test_result1), nrow(peptide_df) * (nrow(peptide_df) + 1) / 2)
  expect_s3_class(test_result1, "data.table")
  expect_s3_class(test_result1, "data.frame")

  test_result2 <- peptide_pairwise_alignment(peptides = peptide_df,
                                             self_comparison = FALSE,
                                             output_str = "tibble")

  expect_equal(nrow(test_result2), nrow(peptide_df) * (nrow(peptide_df) - 1) / 2)
  expect_s3_class(test_result2, "tbl_df")

  expect_equal(unique(test_result1$subject_id), peptide_df$id)
  expect_equal(unique(test_result1$pattern_id), peptide_df$id)
})
#> Test passed
