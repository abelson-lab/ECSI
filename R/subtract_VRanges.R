#' Subtract VRanges objects
#'
#' Identifies features from B that overlap with A, removes the overlapping features from A and returns the remaining portion of A.
#'
#' @param A \code{VRanges} object to filter from and return
#' @param B \code{VRanges} object with alleles to filter out
#' @export
#' @examples
#' \dontrun{
#' variants <- subtract_VRanges(A = variants, B = flagged_alleles)
#' }
#' @return This function returns \code{VRanges} object A with alleles from \code{VRanges} object B filtered out

subtract_VRanges <-
function(A,B) {

  # find matches (specific to chr position strand ref alt)
  m = BiocGenerics::match(A, B)

  # remove matches
  filtered = A[is.na(m)]

  return(filtered)
}
