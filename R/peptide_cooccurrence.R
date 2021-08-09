#' Summarize the co-occurrence of two peptides
#'
#' calculate a contingency table of two peptides' presence in multiple individuals
#' and perform independency test and cooccurrence coefficients
#'
#' @param x a logical vector of peptide 1's presence in multiple samples
#' @param y a logical vector of peptide 2's presence in multiple samples
#' @param coefficient the cooccurrence coefficient to compute.
#'     Options: "phi" (phi correlation coefficient), "prop" (mean proportion of
#'     both peptide being enriched given one of them is enriched), "both"
#'     (both "phi" and "prop").
#'     Default: "phi"
#' @param perform_test a logical indicator of whether statistical tests
#'     should be performed for each comparison.
#'     The test method is specified by \code{test_method}.
#'     Default: \code{FALSE}
#' @param test_method method to perform statistical test.
#'     Options: "fisher" (Fisher's exact test) or "chisq" (Pearson's Chi-squared test).
#'     Default: "fisher"
#'
#' @return if \code{coefficient} = "phi" and \code{perform_test} = \code{FALSE},
#'     return a single number, i.e., the phi coefficient. Otherwise, return
#'     a vector that contains at least these four elements:
#'     TT: frequency of both peptide positive;
#'     TF: frequency of peptide 1 positive and peptide 2 negative;
#'     FT: frequency of peptide 1 genative and peptide 2 position;
#'     FF: frequency of both peptide negative.
#'     Depending on \code{coefficient}, the vector also contains phi coefficient
#'     and/or the mean proportion of cooccurrence. If \code{perform_test} is
#'     \code{TRUE}, the vector also contains the key statistic ("oddsr" - odds ratio,
#'     if using the Fisher's exact test, or "xsq" - chi-squared, if using the Chi-squared
#'     test), and the p-value.
#'
#' @author Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{peptide_pairwise_correlation}},
#'     \code{\link[stats]{fisher.test}},
#'     \code{\link[stats]{chisq.test}}
#' @keywords peptide co-occurrence
#'
#' @examples
#' data(peptide_z)
#' peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)
#' peptide_cooccurrence(x = peptide_bi[, 1], y = peptide_bi[, 2],
#'                      coefficient = "both",
#'                      perform_test = TRUE,
#'                      test_method = "fisher")
#' peptide_cooccurrence(x = peptide_bi[, 1], y = peptide_bi[, 2],
#'                      coefficient = "phi",
#'                      perform_test = TRUE,
#'                      test_method = "chisq")
#' peptide_cooccurrence(x = peptide_bi[, 1], y = peptide_bi[, 2],
#'                      coefficient = "prop",
#'                      perform_test = FALSE)
#' peptide_cooccurrence(x = peptide_bi[, 1], y = peptide_bi[, 2],
#'                      coefficient = "phi",
#'                      perform_test = FALSE)
#'
#' @export
#'
peptide_cooccurrence <- function(x, y,
                                 coefficient = "phi",
                                 perform_test = FALSE,
                                 test_method = "fisher"){

  if(coefficient == "phi" & perform_test == FALSE){

    return(suppressWarnings(cor(x, y, method = "pearson")))

  }else{

    tr <- c(TT = sum(x & y),
            TF = sum(x & (!y)),
            FT = sum((!x) & y),
            FF = sum((!x) & (!y)))

    # Add test results if desired:
    if(perform_test == TRUE & test_method == "fisher"){
      suppressWarnings(fisher_t <- fisher.test(matrix(tr, nrow = 2)))
      tr <- c(tr,
              oddsr = ifelse(is.null(fisher_t$estimate), NA, fisher_t$estimate),
              p = fisher_t$p.value)
    }
    if(perform_test == TRUE & test_method == "chisq"){
      suppressWarnings(chi_t <- chisq.test(matrix(tr, nrow = 2)))
      tr <- c(tr,
              xsq = ifelse(is.null(chi_t$statistic), NA, chi_t$statistic),
              p = chi_t$p.value)
    }

    # Add cooccurrence coefficient if desired
    if(coefficient == "phi"){
      tr <- c(tr, phi = suppressWarnings(cor(x, y, method = "pearson")))
    }
    if(coefficient == "prop"){
      TT_NT <- ifelse(test = (tr["TT"] + tr["FT"]) > 0,
                      yes = tr["TT"] / (tr["TT"] + tr["FT"]),
                      no = 0)
      TT_TN <- ifelse(test = (tr["TT"] + tr["TF"]) > 0,
                      yes = tr["TT"] / (tr["TT"] + tr["TF"]),
                      no = 0)
      tr <- c(tr, mean_prop = mean(TT_NT, TT_TN))
    }
    if(coefficient == "both"){
      TT_NT <- ifelse(test = (tr["TT"] + tr["FT"]) > 0,
                      yes = tr["TT"] / (tr["TT"] + tr["FT"]),
                      no = 0)
      TT_TN <- ifelse(test = (tr["TT"] + tr["TF"]) > 0,
                      yes = tr["TT"] / (tr["TT"] + tr["TF"]),
                      no = 0)
      tr <- c(tr,
              phi = suppressWarnings(cor(x, y, method = "pearson")),
              mean_prop = mean(TT_NT, TT_TN))
    }

    return(tr)
  }
}
