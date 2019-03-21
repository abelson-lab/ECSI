#' Subtract VRanges objects
#'
#' Identifies features from B that overlap with A, removes the overlapping features from A and returns the remaining portion of A. If both A and B are VRanges objects, then this function is allele specific.
#'
#' @param A \code{VRanges} object to filter from and return
#' @param B \code{VRanges} or \code{GRanges} object with alleles to filter out
#' @export
#' @examples
#' \dontrun{
#' variants <- subtract_VRanges(A = variants, B = flagged_alleles)
#' }
#' @return This function returns \code{VRanges} object A with alleles from \code{VRanges} object B filtered out

subtract_VRanges <-
function(A,B) {

  # ensure A is VRanges and B is either VRanges or GRanges
  if(class(A) != "VRanges"){ stop('A must be a VRanges object') }
  if(!(class(B) %in% c("GRanges", "VRanges"))){ stop('B must be either a GRanges or VRanges object') }

  # output if B is VRanges object
  if(class(B) == "VRanges"){
    # find matches (specific to chr position strand ref alt)
    m = BiocGenerics::match(A, B)
    # remove matches
    filtered = A[is.na(m)]

  } else if (class(B) == "GRanges"){
    # remove matches in GRanges
    filtered = A[-S4Vectors::queryHits(GenomicRanges::findOverlaps(A, B))]
  }

  return(filtered)
}
