#' Intersect VRanges objects
#'
#' Identifies features from B that overlap with A, removes the overlapping features from A and returns the remaining portion of A. If both A and B are VRanges objects, then this function is allele specific.
#'
#' @param A \code{VRanges} object to filter from and return
#' @param B \code{VRanges} or \code{GRanges} object with alleles to keep
#' @export
#' @examples
#' \dontrun{
#' variants <- intersect_VRanges(A = variants, B = positions_of_interest)
#' }
#' @return This function returns \code{VRanges} object A retaining only the alleles that overlap with \code{VRanges} object B.

intersect_VRanges <-
function(A,B) {

  # ensure A is VRanges and B is either VRanges or GRanges
  if(class(A) != "VRanges"){ stop('A must be a VRanges object') }
  if(!(class(B) %in% c("GRanges", "VRanges"))){ stop('B must be either a GRanges or VRanges object') }

  # output if B is VRanges object
  if(class(B) == "VRanges"){
    # find matches (specific to chr position strand ref alt)
    m = BiocGenerics::match(A, B)
    # keep matches
    filtered = A[!is.na(m)]

  } else if (class(B) == "GRanges"){
    # keep unique matches in GRanges
    filtered = unique(A[S4Vectors::queryHits(GenomicRanges::findOverlaps(A, B))])

  }

  return(filtered)
}
