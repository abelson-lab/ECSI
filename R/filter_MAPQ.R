#' Filter MAPQ
#'
#' Filter variants to keep those that are equal to or greater than a MAPQ cutoff.
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param MAPQ_cutoff Minimum acceptable MAPQ score for variants to keep
#' @import VariantAnnotation
#' @examples
#' variants <- filter_MAPQ(varscan_output = variants, MAPQ_cutoff = 59)
#' @return This function returns a \code{VRanges} object filtered by MAPQ cutoff

filter_MAPQ <-
function(varscan_output, MAPQ_cutoff = 59) {
  # keep variants above MAQ cutoff
  varscan_output <- varscan_output[which(varscan_output$MapQual2 >= MAPQ_cutoff),]

  return(varscan_output)
}
