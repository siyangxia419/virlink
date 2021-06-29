#' Generate a kmer data frame
#'
#' convert a peptide data frame to kmer data frame
#'
#' @param peptides a tibble of peptides with two columns named "id" and "pep_aa"
#' @param kmer_length length of desired kmer
#' @param overlap number of overlapped residuals between two consecutive kmers
#'
#' @return a data frame of kmers with the desired structure for alignment
#'
#' @author Jennifer L. Remmel, Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{get_kmers}}
#' @keywords peptide sequence manipulation
#'
#' @examples
#' data(peptide_df)
#' kmer_df <- make_kmer_dataframe(peptides = peptide_df)
#' head(kmer_df)
#'
#' @export
#' @import dplyr
#' @import tidyr
#'
make_kmer_dataframe <- function(peptides, kmer_length = 12, overlap = 6){

  ## create a kmer data frame
  dkmers <- peptides %>%
    dplyr::arrange(id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(kmer = list(get_kmers(pep_aa, kmer_length, overlap))) %>%
    tidyr::unnest(cols = kmer) %>% # unnest takes list-cols and makes a row for each element in the list
    dplyr::rename(peptide_id = id, peptide_seq = pep_aa) %>%
    dplyr::group_by(peptide_id) %>%
    dplyr::mutate(id = paste(peptide_id, 1:n(), sep = "_"),
                  start_pos = (1:n() - 1) * (kmer_length - overlap) + 1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(pep_aa = kmer) %>%
    dplyr::relocate(id, pep_aa, start_pos, peptide_id, peptide_seq)

  return(dkmers)
}
