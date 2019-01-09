#' Sequence Context
#'
#' Get sequence context for a \code{VRanges} object of the varscan pileup2cns output
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param context Number of nucleotides to return. For example, 3 would return three nucleotides with the position of the variant being the middle nucleotide.
#' @import VariantAnnotation
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @export
#' @examples
#' variants <- load_as_VRanges(sample_name = "pt123", sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' variants <- sequence_context(varscan_output = variants, context = 3)
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{context}.

### TO DO LIST:
#	- Check that input is VRanges
#	- Accommodate hg38 as well

sequence_context <-
function(varscan_output, context = 3) {

  # Get the substitution and sequence context
  varscan_output$context <- getSeq(BSgenome.Hsapiens.UCSC.hg19, varscan_output + (context - 1)/2)

  # Creating FlankingSeqGroup
      # C > T with 5' A and 3' G would be represented as "A[CT]G".
      # 192 groups in total when using trinucleotide context (recommended)
  varscan_output$FlankingSeqGroup <- paste0(as.character(subseq(varscan_output$context, start = 1, end = (context-1)/2)),
                     "[", ref(varscan_output), alt(varscan_output), "]", as.character(subseq(varscan_output$context, start = (context-1)/2 + 2, end = context)))

  # Remove context to save memory
  varscan_output$context <- NULL

  return(varscan_output)
}
