% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peptide_pairwise_alignment.R
\name{peptide_pairwise_alignment}
\alias{peptide_pairwise_alignment}
\title{Peptide pairwise alignment - vectorized version}
\usage{
peptide_pairwise_alignment(
  peptides,
  id_col = "id",
  seq_col = "pep_aa",
  sub_matrix = "BLOSUM62",
  gap_opening = 10,
  gap_extension = 4,
  align_type = "local",
  self_comparison = TRUE,
  full_align = FALSE,
  other_info = TRUE,
  parallel_ncore = NULL,
  output_str = "data.table"
)
}
\arguments{
\item{peptides}{a tibble of peptide metadata with at least two columns:
peptide id and peptide amino acid sequences.}

\item{id_col}{name of the column that contains peptide id.
Default: "id"}

\item{seq_col}{name of the column that contains peptide amino acid sequences.
Default: "pep_aa"}

\item{sub_matrix}{substitution matrix to use for the alignment.
Default = "BLOSUM62"}

\item{gap_opening}{penalty for starting a gap,
pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
Default = 10}

\item{gap_extension}{penalty for extending a gap,
pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
Default = 4}

\item{align_type}{alignment type,
pass to function \code{\link[Biostrings]{pairwiseAlignment}}.
Options: "global", "local", "overlap", "global-local", and "local-global".
Default: "local"}

\item{self_comparison}{whether alignment with self is performed.
Default = TRUE}

\item{full_align}{whether to keep the full alignment results as a column in the final results.
Default = FALSE}

\item{other_info}{whether to keep peptide information in the output.
Default = TRUE}

\item{parallel_ncore}{number of cores used for internal parallel in \code{data.table}.
Default = NULL}

\item{output_str}{format of the output to return.
Options: "data.table", "data.frame", and "tibble".
Default = "data.table"}
}
\value{
a data frame or data table or tibble with sequence alignment results from all unique pairs of peptides
}
\description{
perform pairwise sequence alignment for a dataset of peptides taking advantange of data.table and vectorized
pairwiseAlignment function in the package "Biostrings"
}
\examples{
data(peptide_df)
pairwise_align <- peptide_pairwise_alignment(peptides = peptide_df)
head(pairwise_align)

}
\references{
pairwiseAlignment: \url{https://bioconductor.org/packages/devel/bioc/vignettes/Biostrings/inst/doc/PairwiseAlignments.pdf}
}
\seealso{
\code{\link[Biostrings]{pairwiseAlignment}}
}
\author{
Jennifer L. Remmel, Siyang Xia \email{sxia@hsph.harvard.edu}
}
\keyword{alignment}
\keyword{peptide}
