#' Summarize the co-occurrence of two peptides
#'
#' calculate a contingency table of two peptides' presence in multiple individuals and
#' perform independency test and phi correlation coefficient
#'
#' @param x1 a logical vector of peptide 1's presence in multiple samples
#' @param x2 a logical vector of peptide 2's presence in multiple samples
#' @param test_method method to perform statistical test, options: "fisher" or "chisq" or "none"
#' @param correlation whether phi correlation is calculated
#'
#' @return a vector that contains 4-7 elements (depending on whether to perform test and correlation):
#'   TT: frequency of both peptide positive
#'   TF: frequency of peptide 1 positive and peptide 2 negative
#'   FT: frequency of peptide 1 genative and peptide 2 position
#'   FF: frequency of both peptide negative
#'   stat: key statistics from the selected test (odds ratio for Fisher's exact test and chi-squared for chi-squared test)
#'   p: p value from the selected test
#'   phi: phi correlation coefficient (calculated by Pearson correlation)
#'
#' @author Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{peptide_pairwise_coocurrence}}
#' @keywords peptide co-occurrence
#'
#' @examples
#' data(peptide_z)
#' peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)
#' cooccur <- peptide_cooccurrence(x1 = peptide_bi[, 1], x2 = peptide_bi[, 2])
#' print(cooccur)
#'
#' @export
#'
peptide_cooccurrence <- function(x1, x2,
                                test_method = "fisher",
                                correlation = TRUE){

  v <- c(TT = sum(x1 & x2),
         TF = sum(x1 & (!x2)),
         FT = sum((!x1) & x2),
         FF = sum((!x1) & (!x2)))

  # Test of association
  if(test_method == "chisq"){         # Chi-square test
    suppressWarnings(chi_t <- chisq.test(matrix(v, nrow = 2)))
    v <- c(v, chi_t$statistic, p = chi_t$p.value)
  }else if(test_method == "fisher"){  # Fisher's exact test
    suppressWarnings(fisher_t <- fisher.test(matrix(v, nrow = 2)))
    v <- c(v, fisher_t$statistic, p = fisher_t$p.value)
  }

  # Phi coefficient: calculated by Pearson correlation, equals to sqrt(chi-square / n)
  if(correlation) v <- c(v, phi = suppressWarnings(cor(x1, x2, method = "pearson")))

  return(v)
}
