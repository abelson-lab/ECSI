#' Sequence Context
#'
#' Get sequence context for a \code{VRanges} object of the varscan pileup2cns output
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param genome Genome of input object, accepts either hg19 or hg38
#' @param context Number of nucleotides to return. For example, 3 would return three nucleotides with the position of the variant being the middle nucleotide.
#' @export
#' @examples
#' \dontrun{
#' variants <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' variants <- sequence_context(varscan_output = variants, context = 3)
#' }
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{context}.

### TO DO LIST:
#	- Check that input is VRanges
#	- Perhaps implement genome auto-detection

sequence_context <-
function(varscan_output, genome = c("hg19", "hg38"), context = 3) {

  if(genome == "hg19"){
    # Get the substitution and sequence context
    varscan_output$context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, varscan_output + (context - 1)/2)
  } else if(genome == "hg38"){
    # Get the substitution and sequence context
    varscan_output$context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, varscan_output + (context - 1)/2)
  }

  # Creating FlankingSeqGroup
      # C > T with 5' A and 3' G would be represented as "A[CT]G".
      # 192 groups in total when using trinucleotide context (recommended)
  varscan_output$FlankingSeqGroup <- paste0(as.character(Biostrings::subseq(varscan_output$context, start = 1, end = (context-1)/2)),
                     "[", VariantAnnotation::ref(varscan_output), VariantAnnotation::alt(varscan_output), "]",
                     as.character(Biostrings::subseq(varscan_output$context, start = (context-1)/2 + 2, end = context)))

  # Remove context to save memory
  varscan_output$context <- NULL

  return(varscan_output)
}
