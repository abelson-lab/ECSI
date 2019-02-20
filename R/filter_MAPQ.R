#' Filter MAPQ
#'
#' Filter variants to keep those that are equal to or greater than a MAPQ cutoff.
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param MAPQ_cutoff_ref Minimum acceptable MAPQ score for reference allele
#' @param MAPQ_cutoff_alt Minimum acceptable MAPQ score for alternate allele
#' @export
#' @examples
#' \dontrun{
#' variants <- filter_MAPQ(varscan_output = variants, MAPQ_cutoff = 59)
#' }
#' @return This function returns a \code{VRanges} object filtered by MAPQ cutoff

filter_MAPQ <-
function(varscan_output, MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59) {

  # remove variant MAPQs below cutoff
  varscan_output <- varscan_output[BiocGenerics::which(varscan_output$MapQual2 >= MAPQ_cutoff_alt),]

  # remove ref MAPQs below cutoff but keep MAPQ==0 (homozygous var)
  varscan_output <- varscan_output[BiocGenerics::which(varscan_output$MapQual1 >= MAPQ_cutoff_ref | varscan_output$MapQual1 == 0),]

  return(varscan_output)
}
