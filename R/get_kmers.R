#' Generate kmers from long peptide sequence
#'
#' given a sequence, return all kmers tiled along the sequence
#'
#' @param sequence sequence of a single peptide
#' @param kmer_length length of desired kmer
#' @param overlap number of overlapped residuals between two consecutive kmers
#'
#' @return a character vector that contains all kmers generated from the input peptide sequence
#'
#' @author Jennifer L. Remmel, Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @seealso \code{\link{make_kmer_dataframe}}
#' @keywords peptide sequence manipulation
#'
#' @examples
#' km <- get_kmers(sequence = "MRSLLFVVGAWVAALVTNLTPDAALASGTTTTAAAGNTSATASPGDNATSIDAGST",
#'                 kmer_length = 12,
#'                 overlap = 6)
#' print(km)
#'
#' @export
#' @import dplyr
#'
get_kmers <- function(sequence, kmer_length = 12, overlap = 6){

  # determine start indices
  starts <- seq(from = 1, to = nchar(sequence), by = (kmer_length - overlap))

  # determine end indices
  ends <- starts + kmer_length - 1

  # slice sequence by start/end indices
  scratch <- tibble::tibble(starts = starts, ends = ends) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(kmers = substr(x = sequence, start = starts, stop = ends)) %>%
    dplyr::filter(nchar(kmers) >= overlap) # remove sequences shorter than the overlap length

  return(scratch$kmers)
}
