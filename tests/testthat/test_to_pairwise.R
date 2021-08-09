library(virlink)

data(peptide_z)

test_that("input d is a numeric data frame and f requires two vectors as inputs", {
  expect_error(to_pairwise(d = cbind.data.frame(rep("a", 10), rep("b", 10)), f = cor))
  expect_warning(to_pairwise(d = matrix(1:50, nrow = 5, ncol = 10), f = cor))
  expect_error(to_pairwise(d = peptide_z, f = dim))
})
#> Test passed

test_that("parallel computing setting", {
  expect_warning(to_pairwise(d = peptide_z, f = cor, mc = -1))
  if(Sys.info()[['sysname']] == "Windows"){
    expect_warning(to_pairwise(d = peptide_z, f = cor, mc = 2))
  }else{
    expect_warning(to_pairwise(d = peptide_z, f = cor, mc = 1000))
  }
})
#> Test passed

test_that("dimension of the results and the names of the two id columns", {
  expect_equal(nrow(to_pairwise(d = peptide_z, f = cor, 
                             unique_pair = TRUE, same_comparison = FALSE)),
               ncol(peptide_z) * (ncol(peptide_z) - 1) / 2)
  
  expect_equal(nrow(to_pairwise(d = peptide_z, f = cor, 
                                unique_pair = TRUE, same_comparison = TRUE)),
               ncol(peptide_z) * (ncol(peptide_z) + 1) / 2)
  
  expect_equal(nrow(to_pairwise(d = peptide_z, f = cor, 
                                unique_pair = FALSE, same_comparison = FALSE)),
               ncol(peptide_z) * (ncol(peptide_z) - 1))
  
  expect_equal(nrow(to_pairwise(d = peptide_z, f = cor, 
                                unique_pair = FALSE, same_comparison = TRUE)),
               ncol(peptide_z) * ncol(peptide_z))
  
  test_result <- to_pairwise(d = peptide_z, f = cor, 
                             unique_pair = TRUE, same_comparison = FALSE)
  expect_equal(colnames(test_result)[1:2], c("id1", "id2"))
  expect_equal(sort(unique(c(test_result$id1, test_result$id2))), 
               sort(colnames(peptide_z)))
  expect_s3_class(test_result, class = "data.frame")
})
#> Test passed
