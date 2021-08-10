#' Peptide pairwise correlation or co-occurrence
#'
#' perform correlation or co-occurrence analysis across all pairs of peptides.
#'     Pairwise comparison were implemented by \code{\link{to_pairwise}}.
#'
#' @param d a data frame or tibble containing z-scores or binary indicator of
#'     peptide enrichment. Columns represent peptides and rows represent samples.
#'     Column names are peptide ids and row names are sample ids.
#' @param analysis_type type of analysis to perform.
#'     Options: "correlation" or "cooccurrence".
#'     If select "correlation", correlation will be calculated on the z-scores, and
#'     \code{cor_method} will be used.
#'     If select "cooccurrence", z-score will be converted to a binary indicator
#'     of enrichment (1/0), and \code{occ_method} and \code{occ_test} will be used.
#'     Default: \code{correlation}
#' @param perform_test a logical indicator of whether statistical tests
#'     should be performed for each comparison.
#'     The test method is specified by \code{cor_method} and \code{occ_method},
#'     depending on \code{analysis_type}.
#'     Default: \code{FALSE}
#' @param cor_method correlation coefficient and test to be computed.
#'     Pass to \code{\link[stats]{cor}} and \code{\link[stats]{cor.test}}
#'     Options: "pearson", "spearman", "kendall".
#'     Default: "pearson"
#' @param occ_method cooccurrence coefficient to be computed.
#'     Pass to \code{\link{peptide_cooccurrence}}.
#'     Options: "phi", "prop", "both".
#'     Default: "phi"
#' @param occ_test cooccurrence tests to be computed.
#'     Pass to \code{\link{peptide_cooccurrence}}.
#'     Options: "fisher" (Fisher's exact test) or "chisq" (Pearson's Chi-squared test).
#'     Default: "fisher"
#' @param p_adjust_method method to adjust p values for multiple comparison.
#'     Pass to \code{\link[stats]{p.adjust}}.
#'     options: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".
#'     Default: "holm"
#' @param z_threshold the z-score threshold above which the epitope is considered as enriched.
#'     Default: 5
#' @param temporal_samples a logical indicator of whether there are multiple
#'     temporal samples from each individual. If \code{TRUE}, an additional data frame
#'     \code{si} is needed to indicate the attribution of samples to individuals.
#'     To account for the temporal structure, z-score difference between each two consecutive
#'     temporal points will be calculated and used for the correlation analysis.
#'     Although the function tolerates cooccurrence analysis with temporal samples,
#'     it is not preferred as changes of binary status are less informative.
#'     To perform cooccurrence analysis, only the change from 0 to 1 between two time points
#'     is considered as 1 in the difference data.
#'     Default: \code{FALSE}
#' @param si a data frame or tibble indicating the temporal structure of the samples.
#'     The first column of \code{si} is the sample id, which should match the row names of \code{d}.
#'     The second column of \code{si} is the individual id.
#'     All other columns will be ignored.
#' @param output_str format of the output to return.
#'     Options: "data.table", "data.frame", and "tibble".
#'     Default: "data.table"
#' @param mc number of cores to use for multi-thread parallel analysis. Pass to
#'     \code{\link[parallel]{mclapply}}.
#'     Deafult: 1
#'
#' @return a data frame/data table/tibble that contains correlation or co-occurence
#'     analysis results of all pairs of peptides
#'
#' @author Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{peptide_coocurrence}}, \code{\link{to_pairwise}}
#' @keywords peptide correlation cooccurrence
#'
#' @examples
#' data(peptide_z)
#' sample_info <- data.frame(id = row.names(peptide_z),
#'                           individual = rep(c("ind1", "ind2"), each = 6))
#'
#' ### 1. correlation analysis with default setting
#' output1 <- peptide_pairwise_correlation(d = peptide_z)
#' head(output1)
#'
#' ### 2. correlation with statistical tests
#' output2 <- peptide_pairwise_correlation(d = peptide_z,
#'                                         si = sample_info,
#'                                         analysis_type = "correlation",
#'                                         perform_test = TRUE,
#'                                         cor_method = "pearson")
#' head(output2)
#'
#' ### 3. correlation with temporal samples
#' output3 <- peptide_pairwise_correlation(d = peptide_z,
#'                                         si = sample_info,
#'                                         analysis_type = "correlation",
#'                                         perform_test = TRUE,
#'                                         cor_method = "pearson",
#'                                         temporal_samples = TRUE)
#' head(output3)
#'
#' ### 4. cooccurrence analysis with both coefficients
#' output4 <- peptide_pairwise_correlation(d = peptide_z,
#'                                         si = sample_info,
#'                                         analysis_type = "cooccurrence",
#'                                         z_threshold = 5,
#'                                         perform_test = FALSE,
#'                                         occ_method = "both",
#'                                         temporal_samples = FALSE)
#' head(output4)
#'
#' ### 5. cooccurrence analysis with both coefficients and fisher's exact test
#' output5 <- peptide_pairwise_correlation(d = peptide_z,
#'                                         si = sample_info,
#'                                         analysis_type = "cooccurrence",
#'                                         z_threshold = 10,
#'                                         perform_test = TRUE,
#'                                         occ_method = "both",
#'                                         occ_test = "fisher",
#'                                         temporal_samples = FALSE)
#' head(output5)
#'
#' ### 6. cooccurrence analysis with only phi coefficient and return as a tibble
#' output6 <- peptide_pairwise_correlation(d = peptide_z,
#'                                         si = sample_info,
#'                                         analysis_type = "cooccurrence",
#'                                         z_threshold = 5,
#'                                         perform_test = FALSE,
#'                                         occ_method = "phi",
#'                                         temporal_samples = FALSE,
#'                                         output_str = "tibble")
#' head(output6)
#'
#' @export
#' @import data.table
#'
peptide_pairwise_correlation <- function(d,
                                         analysis_type    = "correlation", # option: "correlation", "cooccurrence"
                                         perform_test     = FALSE,
                                         cor_method       = "pearson",     # option: "pearson", "spearman"
                                         occ_method       = "phi",         # option: "phi", "prop", "both"
                                         occ_test         = "fisher",      # option: "fisher", "chisq"
                                         p_adjust_method  = "holm",
                                         z_threshold      = 5,
                                         temporal_samples = FALSE,
                                         si               = NULL,
                                         output_str       = "data.table",  # options: "data.table", "data.frame", "tibble"
                                         mc               = 1
){


  # 1. basic messages about the analysis to perform -------------------------

  if(analysis_type == "correlation"){
    message("Correlation of z-scores between each pair of epitopes.")

    if(cor_method == "pearson"){
      message("Correlation method: Pearson")
    }else if(cor_method == "spearman"){
      message("Correlation method: Spearman")
    }else if(cor_method == "kendall"){
      message("Correlation method: Kendall")
    }else{
      stop('Invalid correlation method. Please select from "pearson", "spearman", "kendall".')
    }

    if(perform_test) message("Perform statistical test")

  }else if(analysis_type == "cooccurrence"){
    message("Co-occurrence (binary) analysis between each pair of epitopes.")

    if(occ_method == "phi"){
      message("Cooccurrence coefficient: phi correlation")
    }else if(occ_method == "prop"){
      message("Cooccurrence coefficient: mean proportion of co-occurrence given one peptide is present")
    }else if(occ_method == "both"){
      message("Cooccurrence coefficient: both phi and mean proportion")
    }else{
      stop('Invalid cooccurrence coefficient. Please select from "phi", "prop", "both".')
    }

    if(perform_test){
      if(occ_test == "fisher"){
        message("Perform statistical test: Fisher's exact test")
      }else if(occ_test == "chisq"){
        message("Perform statistical test: Chi-squared test")
      }else{
        stop('Invalid cooccurrence test method Please select from "fisher", "chisq", "both".')
      }
    }

  }else{
    stop('Invalid data type. Please select from "correlation" and "cooccurrence".')
  }

  if(temporal_samples){
    message("Consider multiple temporal samples of each individual")
  }



  # 2. select the function to analyze each pair of epitopes -----------------

  if(analysis_type == "correlation" & perform_test == FALSE){

    fct <- function(v1, v2){
      suppressWarnings(cor(x = v1, y = v2, method = cor_method))
    }

  }else if(analysis_type == "correlation" & perform_test == TRUE){

    fct <- function(v1, v2){
      tr <- suppressWarnings(cor.test(x = v1, y = v2, method = cor_method))
      return(c(tr$estimate, tr$statistic, p = tr$p.value))
    }

  }else if(analysis_type == "cooccurrence"){

    fct <- function(v1, v2) peptide_cooccurrence(x = v1, y = v2,
                                                 coefficient = occ_method,
                                                 perform_test = perform_test,
                                                 test_method = occ_test)
  }



  # 3. pairwise analysis ----------------------------------------------------

  ### convert the z-score data to binary occurrence data if necessary
  if(analysis_type == "cooccurrence"){
    if(!all(unlist(d) %in% c(0, 1))) d <- as.data.frame((d > z_threshold) * 1)
  }


  ### pairwise analysis
  if(temporal_samples == FALSE){  # do not consider temporal samples

    output <- to_pairwise(d = d,
                          f = fct,
                          unique_pair = TRUE,
                          same_comparison = FALSE,
                          mc = mc)

  }else{  # use temporal samples

    # check if "si" contains the correct information
    if(length(intersect(unlist(si[, 1]), row.names(d))) < 2){
      stop("Too few samples are included in the sample information data.")
    }

    # colnames of d and si
    dt_name <- names(d)

    # convert the original data and the sample information (si) to data tables
    d_t <- data.table::setDT(data.table::copy(d), keep.rownames = "id")
    s_i <- data.table::setDT(data.table::copy(si[, 1:2]))
    data.table::setnames(s_i, new = c("id", "individual"))

    # combine them
    d_t <- d_t[s_i, on = "id"]
    data.table::setcolorder(d_t, neworder = c("id", "individual"))

    # calculate changes between every two consecutive time point for each individual and each epitope
    d_t[, (dt_name) := .SD - data.table::shift(.SD, give.names = TRUE), by = .(individual), .SDcols = dt_name]
    d_t <- na.omit(d_t)

    # check if "si" contains the correct information
    if(d_t[, .N] < 2){
      stop("Too few samples left after calculating changes between time points.")
    }else if(d_t[, .N] < 10){
      warning("Fewer than 10 samples left after calculating changes between time points.")
    }

    # remove sample information and keep only the data
    d_t[, c("id", "individual") := NULL]

    # when the original data is presence/absence (binary), calculating difference will produce -1
    if(analysis_type == "cooccurrence"){
      d_t[d_t < 0] <- 0
      warning("Cooccurrence analysis is not prefered for temporal data. Correlation analysis with z-scores is more informative.")
    }

    # pairwise analysis
    output <- to_pairwise(d = d_t,
                          f = fct,
                          unique_pair = TRUE,
                          same_comparison = FALSE,
                          mc = mc)

  }


  ### convert the output to data table
  if(!data.table::is.data.table(output)) output <- data.table::setDT(output)


  ### fix column name
  if("temp" %in% colnames(output)) data.table::setnames(output, old = "temp", new = "cor")


  ### adjust p-values if performed statistical analysis
  if(perform_test){
    output[, p_adjust := p.adjust(p, method = p_adjust_method)]
  }



  # 4. Return the results -------------------------------------------------

  if(output_str == "data.table"){
    return(output)
  }else if(output_str == "data.frame"){
    return(as.data.frame(output))
  }else if(output_str == "tibble"){
    return(tibble::as_tibble(output))
  }else{
    warning('Invalid output structure. Result returned as a "data.table".')
    return(output)
  }

}
