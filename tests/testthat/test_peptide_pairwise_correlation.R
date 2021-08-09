library(virlink)

data(peptide_z)
sample_info <- data.frame(id = row.names(peptide_z),
                          individual = rep(c("ind1", "ind2"), each = 6))

test_that("produce message for each run", {
  expect_message(peptide_pairwise_correlation(d = peptide_z))
})
#> Test passed

test_that("input d is a numeric data frame", {
  expect_error(expect_message(peptide_pairwise_correlation(d = cbind.data.frame(rep("a", 10), rep("b", 10)))))
})
#> Test passed

test_that("parallel computing setting", {
  expect_warning(expect_message(peptide_pairwise_correlation(d = peptide_z, mc = -1)))
  if(Sys.info()[['sysname']] == "Windows"){
    expect_warning(expect_message(peptide_pairwise_correlation(d = peptide_z, mc = 2)))
  }else{
    expect_warning(expect_message(peptide_pairwise_correlation(d = peptide_z, mc = 1000)))
  }
})
#> Test passed

test_that("input options", {
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "none")))
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "correlation", 
                                                           perform_test = TRUE, 
                                                           cor_method = "phi")))
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           occ_method = "cor")))
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           occ_test = "cor")))
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           p_adjust_method = "xx")))
  expect_warning(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           output_str = "matrix")))
})
#> Test passed

test_that("dimension of the results and the names of the two id columns", {
  expect_equal(nrow(expect_message(peptide_pairwise_correlation(d = peptide_z))),
               ncol(peptide_z) * (ncol(peptide_z) - 1) / 2)
  
  test_result <- expect_message(peptide_pairwise_correlation(d = peptide_z))
  expect_equal(colnames(test_result)[1:2], c("id1", "id2"))
  expect_equal(sort(unique(c(test_result$id1, test_result$id2))), 
               sort(colnames(peptide_z)))
  expect_s3_class(test_result, class = "data.frame")
})
#> Test passed

test_that("handling temporal samples", {
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "correlation", 
                                                           perform_test = TRUE, 
                                                           temporal_samples = TRUE)))
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "correlation", 
                                                           perform_test = TRUE, 
                                                           temporal_samples = TRUE,
                                                           si = sample_info[, 2])))
  expect_warning(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                             analysis_type = "cooccurrence", 
                                                             perform_test = TRUE, 
                                                             temporal_samples = TRUE, 
                                                             si = sample_info)))
  
  sample_info_wrong <- sample_info
  sample_info_wrong$id <- paste0(sample_info_wrong$id, "x")
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           temporal_samples = TRUE,
                                                           si = sample_info_wrong)))
  
  sample_info_wrong <- sample_info
  sample_info_wrong$individual <- sample_info_wrong$id
  expect_error(expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                           analysis_type = "cooccurrence", 
                                                           perform_test = TRUE, 
                                                           temporal_samples = TRUE,
                                                           si = sample_info_wrong)))
  
  expect_false(identical(x = expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                                         analysis_type = "correlation", 
                                                                         perform_test = TRUE, 
                                                                         temporal_samples = TRUE, 
                                                                         si = sample_info)),
                         y = expect_message(peptide_pairwise_correlation(d = peptide_z, 
                                                                         analysis_type = "correlation", 
                                                                         perform_test = TRUE, 
                                                                         temporal_samples = FALSE))))
})
#> Test passed
