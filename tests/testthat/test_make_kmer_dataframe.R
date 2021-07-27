library(virlink)

data(peptide_df)
test_kmer_length <- 9
test_overlap <- 4

test_that("input is a data frame with the first two columns being 'id' and 'pep_aa'", {
  expect_error(make_kmer_dataframe(peptides = peptide_df[, -1]))
  expect_error(make_kmer_dataframe(peptides = peptide_df[, -2]))
  expect_error(make_kmer_dataframe(peptides = "NNNNNNN"))
})
#> Test passed
test_that("output has the right dimensions", {
  test_result <- make_kmer_dataframe(peptides = peptide_df,
                                     kmer_length = test_kmer_length,
                                     overlap = test_overlap)

  expect_equal(ncol(test_result) - ncol(peptide_df), 3)
  expect_equal(sum(ceiling((nchar(peptide_df$pep_aa) - test_overlap) / (test_kmer_length - test_overlap))),
               nrow(test_result))
  expect_equal(unique(test_result$peptide_seq), peptide_df$pep_aa)
})
#> Test passed
