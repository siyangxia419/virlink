library(virlink)

test_aa <- "MRSLLFVVGAWVAALVTNLTPDAALASGTTTTAAAGNTSATASPGDNATSIDAGST"
test_kmer_length <- 8
test_overlap <- 4

test_that("input is a single character string", {
  expect_error(get_kmers(sequence = c(test_aa, test_aa),
                         kmer_length = test_kmer_length,
                         overlap = test_overlap))
  expect_error(get_kmers(sequence = NA,
                         kmer_length = test_kmer_length,
                         overlap = test_overlap))
})
#> Test passed
test_that("returns a character vector of expected length", {
  test_result <- get_kmers(sequence = test_aa,
                           kmer_length = test_kmer_length,
                           overlap = test_overlap)

  expect_equal(max(nchar(test_result)), test_kmer_length)
  expect_true(min(nchar(test_result)) >= test_overlap)
  expect_equal(sum(nchar(test_result) - test_overlap) + test_overlap, nchar(test_aa))
})
#> Test passed 🎉
