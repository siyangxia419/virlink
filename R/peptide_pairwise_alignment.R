#' Peptide pairwise alignment - vectorized version
#'
#' perform pairwise sequence alignment for a dataset of peptides taking advantange of data.table and vectorized
#' pairwiseAlignment function in the package "Biostrings"
#'
#' @param peptides a tibble of peptide metadata with at least two columns:
#'     peptide id and peptide amino acid sequences.
#' @param id_col name of the column that contains peptide id.
#'     Default: "id"
#' @param seq_col name of the column that contains peptide amino acid sequences.
#'     Default: "pep_aa"
#' @param sub_matrix substitution matrix to use for the alignment.
#'     Default = "BLOSUM62"
#' @param gap_opening penalty for starting a gap,
#'     pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
#'     Default = 10
#' @param gap_extension penalty for extending a gap,
#'     pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
#'     Default = 4
#' @param align_type alignment type,
#'     pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
#'     Options: "global", "local", "overlap", "global-local", and "local-global".
#'     Default: "local"
#' @param self_comparison whether alignment with self is performed.
#'     Default = TRUE
#' @param full_align whether to keep the full alignment results as a column in the final results.
#'     Default = FALSE
#' @param other_info whether to keep peptide information in the output.
#'     Default = TRUE
#' @param parallel_ncore number of cores used for internal parallel in \code{data.table}.
#'     Default = NULL
#' @param output_str format of the output to return.
#'     Options: "data.table", "data.frame", and "tibble".
#'     Default = "data.table"
#'
#' @return a data frame or data table or tibble with sequence alignment results from all unique pairs of peptides
#'
#' @author Jennifer L. Remmel, Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link[Biostrings]{pairwiseAlignment}}
#' @references pairwiseAlignment: \url{https://bioconductor.org/packages/devel/bioc/vignettes/Biostrings/inst/doc/PairwiseAlignments.pdf}
#' @keywords peptide alignment
#'
#' @examples
#' data(peptide_df)
#' pairwise_align <- peptide_pairwise_alignment(peptides = peptide_df)
#' head(pairwise_align)
#'
#' @export
#' @import data.table
#' @import tibble
#'
peptide_pairwise_alignment <- function(peptides,
                                       id_col          = "id",
                                       seq_col         = "pap_aa",
                                       sub_matrix      = "BLOSUM62",
                                       gap_opening     = 10,
                                       gap_extension   = 4,
                                       align_type      = "local",
                                       self_comparison = TRUE,
                                       full_align      = FALSE,
                                       other_info      = TRUE,
                                       parallel_ncore  = NULL,
                                       output_str      = "data.table"){


  # 1. Set the number of cores used by data.table if specified by us --------

  if(!is.null(parallel_ncore)){
    cl <- parallel::detectCores()  # number of clusters available
    if(parallel_ncore > (cl - 1)){
      parallel_ncore <- cl - 1
      print("Required number of cores exceed available cores. Too greedy!")
    }
    data.table::setDTthreads(threads = parallel_ncore)
  }



  # 2. Prepare the data.table -----------------------------------------------

  # convert the data frame to data table
  p <- data.table::setDT(data.table::copy(peptides))
  data.table::setnames(p, old = id_col, new = "id")
  data.table::setnames(p, old = seq_col, new = "pep_aa")
  data.table::setkey(p, id)

  # columns to keep in the final results
  if(other_info){
    keeped_name <- names(p)
  }else{
    keeped_name <- c("id", "pep_aa")
  }

  # copy the dataset for pattern and subject
  p1 <- data.table::copy(p[, keeped_name, with = FALSE])
  data.table::setnames(p1, names(p1), paste0("subject_", names(p1)))
  p2 <- data.table::copy(p[, keeped_name, with = FALSE])
  data.table::setnames(p2, names(p2), paste0("pattern_", names(p2)))

  # convert peptide sequence to AAString
  p[, pep_aas := sapply(pep_aa, Biostrings::AAString)]

  # add the id thresholds for the pairwiseAlignment_wrap2 function
  p[, id_cut := id]



  # 3. Pairwise alignment ---------------------------------------------------

  # wrapper function for pariwiseAlignment
  if(full_align){
    pairwiseAlignment_wrap2 <- function(subject_pep_aas, id_threshold){
      alignment <- Biostrings::pairwiseAlignment(pattern = p[id >= id_threshold, pep_aas],
                                                 subject = subject_pep_aas,
                                                 type = align_type,
                                                 gapOpening = gap_opening,
                                                 gapExtension = gap_extension,
                                                 substitutionMatrix = sub_matrix)
      alignment_list <- sapply(alignment, FUN = function(x) list(x))
      score          <- Biostrings::score(alignment)
      string_compare <- Biostrings::compareStrings(alignment)
      nchar          <- Biostrings::nchar(alignment)
      matches        <- Biostrings::nmatch(alignment)
      mismatches     <- Biostrings::nmismatch(alignment)
      return(list(pattern_id     = p[id >= id_threshold, id],
                  alignment      = alignment_list,
                  string_compare = string_compare,
                  score          = score,
                  nchar          = nchar,
                  matches        = matches,
                  mismatches     = mismatches))
    }
  }else{
    pairwiseAlignment_wrap2 <- function(subject_pep_aas, id_threshold){
      alignment <- Biostrings::pairwiseAlignment(pattern = p[id >= id_threshold, pep_aas],
                                                 subject = subject_pep_aas,
                                                 type = align_type,
                                                 gapOpening = gap_opening,
                                                 gapExtension = gap_extension,
                                                 substitutionMatrix = sub_matrix)
      score          <- Biostrings::score(alignment)
      string_compare <- Biostrings::compareStrings(alignment)
      nchar          <- Biostrings::nchar(alignment)
      matches        <- Biostrings::nmatch(alignment)
      mismatches     <- Biostrings::nmismatch(alignment)
      return(list(pattern_id     = p[id >= id_threshold, id],
                  string_compare = string_compare,
                  score          = score,
                  nchar          = nchar,
                  matches        = matches,
                  mismatches     = mismatches))
    }
  }

  # pairwise alignment
  palign <- p[, pairwiseAlignment_wrap2(subject_pep_aas = pep_aas, id_threshold = id_cut), by = id]



  # 4. Clean up the alignment results ---------------------------------------

  # remove self alignments if needed
  if(!self_comparison) palign <- palign[id != pattern_id,]

  # add number of gaps
  palign[, ":="(gaps = stringr::str_count(string_compare, "[\\+, \\-]"))]

  # add useful peptide information
  data.table::setnames(palign, old = "id", new = "subject_id")
  palign <- palign[p1, on = "subject_id", nomatch = NULL][p2, on = "pattern_id", nomatch = NULL]

  # set keys
  data.table::setkey(x = palign, subject_id, pattern_id)
  data.table::setcolorder(x = palign,
                          neworder = c("subject_id", "pattern_id",
                                       "subject_pep_aa", "pattern_pep_aa",
                                       "string_compare",
                                       "nchar", "matches", "mismatches", "gaps", "score"))


  # Return the results
  if(output_str == "data.table"){
    return(palign)
  }else if(output_str == "data.frame"){
    return(as.data.frame(palign))
  }else if(output_str == "tibble"){
    return(tibble::as_tibble(palign))
  }else{
    return(palign)
    warning('Invalid output structure. Result returned as a "data.table".')
  }
}
