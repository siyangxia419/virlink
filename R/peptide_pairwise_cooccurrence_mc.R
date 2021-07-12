#' Peptide pairwise co-occurrence
#'
#' perform co-occurrence analysis across all pairs of peptides
#'
#' @param d a data frame or tibble with columns indicating peptides and rows indicating samples.
#'     column names are peptide id and row names are sample id
#' @param test_method method to perform statistical test, options: "fisher" or "chisq" or "none"
#' @param correlation whether phi correlation is calculated
#' @param p_adjust_method method to adjust multiple test p values
#'     options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @param parallel_ncore number of cores used for internal parallel in data.table
#'
#' @return a tibble contains co-occurence analysis results of all pairs of peptides
#'
#' @author Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{peptide_coocurrence}}
#' @keywords peptide co-occurrence
#'
#' @examples
#' data(peptide_z)
#' peptide_bi <- as.data.frame((peptide_z > 3.5) * 1)
#' pairwise_cooccur <- peptide_pairwise_cooccurrence(d = peptide_bi)
#' head(pairwise_cooccur)
#'
#' @export
#' @import tibble
#' @import dplyr
#'
peptide_pairwise_cooccurrence_mc <- function(d,
                                             test_method = "fisher",
                                             correlation = TRUE,
                                             p_adjust_method = "holm",
                                             parallel_ncore = 1){

  # number of cores to use
  max_n_core <- parallel::detectCores()  # number of clusters available
  if(parallel_ncore >= max_n_core){
    parallel_ncore <- max_n_core - 1
    print("Required number of cores exceed available cores. Too greedy!")
  }

  co_occur <- parallel::mclapply(X = d,
                                 mc.cores = parallel_ncore,
                                 FUN = function(x){
                                   sapply(d, function(y) peptide_cooccurrence(x1 = x, x2 = y,
                                                                              test_method = test_method,
                                                                              correlation = correlation)) %>%
                                     t() %>%
                                     as.data.frame() %>%
                                     tibble::rownames_to_column(var = "id2") %>%
                                     tibble::as_tibble()
                                 }) %>%
    dplyr::bind_rows(.id = "id1") %>%
    dplyr::mutate(id1 = readr::parse_integer(id1),
                  id2 = readr::parse_integer(id2)) %>%
    dplyr::filter(id1 < id2) %>%
    dplyr::mutate(method = test_method,
                  p_adj = p.adjust(p, method = p_adjust_method),
                  prop_1_in_2 = TT / (TT + FT),
                  prop_2_in_1 = TT / (TT + TF),
                  mean_prop = (prop_1_in_2 + prop_2_in_1) / 2)

  return(co_occur)
}
