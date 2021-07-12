#' Peptide pairwise alignment - element-wise version
#'
#' perform pairwise sequence alignment for a dataset of peptides using dplyr and vectorized
#' pairwiseAlignment function in the package "Biostrings"
#'
#' @param peptides a tibble of peptides with two columns named "id" and "pep_aa"
#' @param sub_matrix substitution matrix to use for the alignment, default = "BLOSUM50"
#' @param gap_opening penalty for starting a gap in the pairwiseAlignment function, default = 10
#' @param gap_extension penalty for extending a gap in the pairwiseAlignment function, default = 4
#' @param align_type alignment type in the pairwiseAlignment function,
#'     options including "local" (default), "global", "overlap", "global-local", and "local-global".
#' @param self_comparison whether alignment with self is performed, default = TRUE
#' @param full_align whether to keep the full alignment results as a column in the final results, default = FALSE
#' @param other_info whether to keep peptide information in the output, default = TRUE
#' @param parallel_ncore number of cores used for internal parallel in data.table
#'
#' @return a tibble with sequence alignment results from all unique pairs of peptides
#'
#' @author Jennifer L. Remmel, Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{peptide_comparisons_vctz}}
#' @references pairwiseAlignment: \url{https://bioconductor.org/packages/devel/bioc/vignettes/Biostrings/inst/doc/PairwiseAlignments.pdf}
#' @keywords peptide alignment
#'
#' @examples
#' data(peptide_df)
#' pairwise_align1 <- peptide_comparisons_elmt(peptides = peptide_df)
#' pairwise_align2 <- peptide_comparisons_elmt(peptides = peptide_df, parallel_ncore = 4)
#' head(pairwise_align1)
#' head(pairwise_align2)
#'
#' @import tibble
#' @import dplyr
#'
peptide_comparisons_elmt <- function(peptides,
                                     sub_matrix      = "BLOSUM50",
                                     gap_opening     = 10,
                                     gap_extension   = 4,
                                     align_type      = "local",
                                     self_comparison = TRUE,
                                     full_align      = FALSE,
                                     other_info      = TRUE,
                                     parallel_ncore  = 1){


  # 1. Preparation ----------------------------------------------------------

  # arrange peptides by id
  peptides <- peptides %>% dplyr::arrange(id)

  # pattern/subject
  patterns <- peptides %>%
    setNames(paste0("pattern_", names(.))) %>%
    dplyr::rename(pattern = pattern_pep_aa)
  if(!other_info) patterns <- patterns %>% dplyr::select(pattern_id, pattern)

  subjects <- peptides %>%
    setNames(paste0("subject_", names(.))) %>%
    dplyr::rename(subject = subject_pep_aa)
  if(!other_info) subjects <- subjects %>% dplyr::select(subject_id, subject)



  # 2. create a wrapper function for pairwiseAlignment ----------------------

  if(full_align){
    pairwiseAlignment_wrap <- function(pattern, subject){
      alignment <- Biostrings::pairwiseAlignment(pattern = Biostrings::AAString(pattern),
                                                 subject = Biostrings::AAString(subject),
                                                 type = align_type,
                                                 gapOpening = gap_opening,
                                                 gapExtension = gap_extension,
                                                 substitutionMatrix = sub_matrix)
      score          <- Biostrings::score(alignment)
      string_compare <- Biostrings::compareStrings(alignment)
      nchar          <- Biostrings::nchar(alignment)
      matches        <- Biostrings::nmatch(alignment)
      mismatches     <- Biostrings::nmismatch(alignment)
      return(list(alignment      = alignment,
                  string_compare = string_compare,
                  score          = score,
                  nchar          = nchar,
                  matches        = matches,
                  mismatches     = mismatches))
    }
  }else{
    pairwiseAlignment_wrap <- function(pattern, subject){
      alignment <- Biostrings::pairwiseAlignment(pattern = Biostrings::AAString(pattern),
                                                 subject = Biostrings::AAString(subject),
                                                 type = align_type,
                                                 gapOpening = gap_opening,
                                                 gapExtension = gap_extension,
                                                 substitutionMatrix = sub_matrix)
      score          <- Biostrings::score(alignment)
      string_compare <- Biostrings::compareStrings(alignment)
      nchar          <- Biostrings::nchar(alignment)
      matches        <- Biostrings::nmatch(alignment)
      mismatches     <- Biostrings::nmismatch(alignment)
      return(list(string_compare = string_compare,
                  score          = score,
                  nchar          = nchar,
                  matches        = matches,
                  mismatches     = mismatches))
    }
  }



  # 3. Pairwise alignment ---------------------------------------------------

  if(parallel_ncore > 1){

    # prepare for parallel computing
    cl <- parallel::detectCores()  # number of clusters available
    if(parallel_ncore > (cl - 1)){
      parallel_ncore <- cl - 1
      print("Required number of cores exceed available cores. Too greedy!")
    }

    # create clusters
    cluster <- multidplyr::new_cluster(n = parallel_ncore)

    # prepare the data frame
    peptide_pairwise <- peptides %>%
      dplyr::mutate(subject_id = id, subject = pep_aa) %>%
      dplyr::rename(pattern_id = id, pattern = pep_aa) %>%
      tidyr::expand(subject_id, pattern_id) %>%
      dplyr::left_join(subjects, by = "subject_id") %>%
      dplyr::left_join(patterns, by = "pattern_id") %>%
      dplyr::filter(if(self_comparison) pattern_id >= subject_id else pattern_id > subject_id)

    # partition the data frame into roughly even chunks
    peptide_pairwise$group <- rep(1:parallel_ncore, length.out = nrow(peptide_pairwise))
    peptide_pair_partition <- peptide_pairwise %>%
      dplyr::group_by(group) %>%
      multidplyr::partition(cluster = cluster)

    # Register library, function and variable
    multidplyr::cluster_library(cluster,
                                packages = c("Biostrings", "BiocGenerics", "dplyr", "tidyr"))
    multidplyr::cluster_copy(cluster,
                             c("pairwiseAlignment_wrap",
                               "align_type", "gap_opening", "gap_extension", "sub_matrix"))

    # peptide alignment:
    peptide_alignments <- peptide_pair_partition %>%
      dplyr::mutate(alignment = purrr::pmap(.l = list(pattern, subject), .f = pairwiseAlignment_wrap)) %>%
      dplyr::collect() %>%
      dplyr::ungroup() %>%
      dplyr::arrange(subject_id, pattern_id) %>%
      tidyr::unnest_wider(col = alignment) %>%
      dplyr::mutate(gaps = stringr::str_count(string_compare, "[\\+, \\-]")) %>%
      dplyr::rename(subject_pep_aa = subject, pattern_pep_aa = pattern) %>%
      dplyr::relocate(subject_id, pattern_id, subject_pep_aa, pattern_pep_aa,
                      string_compare, nchar, matches, mismatches, gaps, score) %>%
      dplyr::select(-group)

  }else{

    # perform sequential alignment using only one core
    peptide_alignments <- peptides %>%
      dplyr::mutate(subject_id = id, subject = pep_aa) %>%
      dplyr::rename(pattern_id = id, pattern = pep_aa) %>%
      tidyr::expand(pattern_id, subject_id) %>%
      dplyr::left_join(patterns, by = "pattern_id") %>%
      dplyr::left_join(subjects, by = "subject_id") %>%
      dplyr::filter(if(self_comparison) pattern_id >= subject_id else pattern_id > subject_id) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(alignment = purrr::pmap(.l = list(pattern, subject), .f = pairwiseAlignment_wrap)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(subject_id, pattern_id) %>%
      tidyr::unnest_wider(col = alignment) %>%
      dplyr::mutate(gaps = stringr::str_count(string_compare, "[\\+, \\-]")) %>%
      dplyr::rename(subject_pep_aa = subject, pattern_pep_aa = pattern) %>%
      dplyr::relocate(subject_id, pattern_id, subject_pep_aa, pattern_pep_aa,
                      string_compare, nchar, matches, mismatches, gaps, score)
  }

  return(peptide_alignments)
}
