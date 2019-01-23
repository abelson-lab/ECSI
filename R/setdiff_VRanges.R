#' Setdiff with VRanges objects
#'
#' Filter variants to keep those that are equal to or greater than a MAPQ cutoff.
#'
#' @param x \code{VRanges} object to filter from and return
#' @param y \code{VRanges} object with alleles to filter out
#' @examples
#' \dontrun{
#' variants <- setdiff_VRanges(x = variants, y = flagged_alleles)
#' }
#' @return This function returns \code{VRanges} object x with alleles from \code{VRanges} object y filtered out

setdiff_VRanges <-
function(x,y) {

  # find matches (specific to chr position strand ref alt)
  m = BiocGenerics::match(x, y)

  # remove matches
  filtered = x[is.na(m)]

  return(filtered)
}
