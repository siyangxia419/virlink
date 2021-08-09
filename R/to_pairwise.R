#' Enable pairwise analysis with customized function
#'
#' given a data frame, perform pairwise analysis with any function on all pairs of columns
#'
#' @param d a data frame or tibble or data table
#' @param f a function to perform some analysis. The inputs of the function should be two numeric vectors.
#' @param unique_pair a binary indicator of whether only unique pairs should be kept.
#'     If \code{TRUE}, "B-A" pair will be removed if "A-B" pair is present.
#' @param same_comparison a binary indicator of whether to keep self pairs (i.e., "A-A").
#'     If \code{TRUE}, self pairs will be kept.
#' @param mc number of cores to use for multi-thread parallel analysis. Pass to
#'     \code{\link[parallel]{mclapply}}.
#'     Deafult: 1
#'
#' @return a data table contains analysis results of all pairs of columns.
#'     The first two columns ("id1" and "id2") contain the two column names in \code{d}
#'     that were included in each pairwise analysis.
#'
#' @author Siyang Xia \email{sxia@@hsph.harvard.edu}
#' @keywords column pairwise analysis
#'
#' @examples
#' data(peptide_z)
#' peptide_bi <- as.data.frame((peptide_z > 5) * 1)
#' output <- to_pairwise(d = peptide_bi, f = peptide_cooccurrence)
#' head(output)
#'
#' @export
#' @import data.table
#'
to_pairwise <- function(d, f,
                        unique_pair = TRUE,
                        same_comparison = FALSE,
                        mc = 1){

  # # convert the input data to data table:
  # if(!data.table::is.data.table(d)) d <- data.table::setDT(d)

  # check if d is a data frame
  if(!is.data.frame(d)){
    d <- as.data.frame(d)
    warning("The input d is not a data frame. Convert it to data frame.")
  }



  # 0. prepare parallel computing -------------------------------------------

  # check if mc is a integer
  if(!is.integer(mc)) mc <- floor(mc)

  # check if mc is smaller than 0
  if(mc <= 0){
    mc <- 1
    warning("Number of cores specified is smaller than 1. Use 1 core as default.")
  }

  # check if mc is larger than the maximum number of cores available
  max_cores <- parallel::detectCores() - 1 # leave one core out
  if(mc > max_cores){
    mc <- max_cores
    warning("Number of cores specified is larger than the maximum number of cores. Use the maximum number of cores instead.")
  }

  # if the system is a windows, mclapply does not work
  if(Sys.info()[['sysname']] == "Windows" & mc > 1){
    mc <- 1
    warning("Windows system does not support function mclapply. The analysis is done without parallel.")
  }



  # 1. using only a single thread -------------------------------------------

  if(mc == 1){
    # pairwise analysis:
    output <- data.table::rbindlist(
      l = lapply(d, function(x){

        # inner cycle:
        temp <- sapply(d, function(y) f(x, y))

        # format the sapply output as a data.table
        if(is.null(dim(temp))){
          temp <- data.table::setDT(data.frame(temp),
                                    keep.rownames = TRUE)
        }else{
          temp <- data.table::data.table(t(temp),
                                         keep.rownames = TRUE,
                                         stringsAsFactors = FALSE)
        }

      }), use.names = TRUE, idcol = "id1"
    )
  }

  # 2. parallel with multiple threads ---------------------------------------

  # avoid using data.table inside mclapply as both uses multiple threads
  if(mc > 1){

    # pairwise analysis:
    output <- data.table::rbindlist(
      l = parallel::mclapply(d, mc.cores =  mc, FUN = function(x){

        # inner cycle:
        temp <- sapply(d, function(y) f(x, y))

        # format the sapply output as a data.table
        if(is.null(dim(temp))){
          temp <- tibble::rownames_to_column(data.frame(temp), var = "rn")
        }else{
          temp <- tibble::rownames_to_column(as.data.frame(t(temp)), var = "rn")
        }

      }), use.names = TRUE, idcol = "id1"
    )

  }



  # 3. Further formatting ---------------------------------------------------

  # convert the output to data table
  if(!data.table::is.data.table(output)) output <- data.table::setDT(output)

  # if the function f returns a list, the resulting output will have multiple columns
  # as "list". To unlist these columns:
  if(any(sapply(output, is.list))){
    output[, colnames(output) := lapply(.SD, unlist), .SDcols = colnames(output)]
  }

  # fix column names:
  data.table::setnames(output, old = "rn", new = "id2")

  # remove duplicated pairs and self pairs if desired:
  output[, ":="(id1 = ordered(id1, levels = colnames(d)),
                id2 = ordered(id2, levels = colnames(d)))]

  if(unique_pair) output <- output[id1 <= id2]
  if(!same_comparison) output <- output[id1 != id2]

  output[, ":="(id1 = as.character(id1),
                id2 = as.character(id2))]


  return(output)
}
