#' Filter MAPQ
#'
#' Filter variants to keep those that are equal to or greater than a MAPQ cutoff.
#'
#' @param sample \code{VRanges} object of the varscan pileup2cns output
#' @param MAPQ_cutoff_ref Minimum acceptable MAPQ score for reference allele
#' @param MAPQ_cutoff_alt Minimum acceptable MAPQ score for alternate allele
#' @export
#' @examples
#' \dontrun{
#' variants <- filter_MAPQ(sample = variants, MAPQ_cutoff = 59)
#' }
#' @return This function returns a \code{VRanges} object filtered by MAPQ cutoff

filter_MAPQ <-
function(sample, MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59) {

  # Check that input is VRanges
  if(class(sample) != "VRanges"){ stop('Input "sample" needs to be VRanges object') }

  # remove variant MAPQs below cutoff
  sample <- sample[BiocGenerics::which(sample$MapQual2 >= MAPQ_cutoff_alt),]

  # remove ref MAPQs below cutoff but keep MAPQ==0 (homozygous var)
  sample <- sample[BiocGenerics::which(sample$MapQual1 >= MAPQ_cutoff_ref | sample$MapQual1 == 0),]

  return(sample)
}
