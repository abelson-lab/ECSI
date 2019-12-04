#' Sequence Context
#'
#' Get sequence context for a \code{VRanges} object of the varscan pileup2cns output
#'
#' @param sample \code{VRanges} object of the varscan pileup2cns output
#' @param genome Genome of input object, either hg19 or hg38
#' @param context Number of nucleotides to return. For example, 3 would return three nucleotides with the position of the variant being the middle nucleotide.
#' @export
#' @examples
#' \dontrun{
#' variants <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' variants <- sequence_context(sample = variants, context = 3)
#' }
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{context}.

### TO DO LIST:
#	- Check that input is VRanges
#	- Perhaps implement genome auto-detection

sequence_context <-
function(sample, genome = c("hg19", "hg38"), context = 3) {

  # Check that genome was specified by the user
  if(length(genome) > 1){ stop("Need to specify genome, either hg19 or hg38.") }
  # Check that input is VRanges
  if(class(sample) != "VRanges"){ stop('Input "sample" needs to be VRanges object') }

  if(genome == "hg19"){
    # Get the substitution and sequence context
    sample$context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, sample + (context - 1)/2)
  } else if(genome == "hg38"){
    # Get the substitution and sequence context
    sample$context <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, sample + (context - 1)/2)
  } else {
    stop(paste("Genome",genome,"not recognized. Genome must be either hg19 or hg38."))
  }

  # Creating FlankingSeqGroup
  # C > T with 5' A and 3' G would be represented as "A[CT]G".
  # 192 groups in total when using trinucleotide context (recommended)
  sample$FlankingSeqGroup <- paste0(as.character(Biostrings::subseq(sample$context, start = 1, end = (context-1)/2)),
                                    "[", VariantAnnotation::ref(sample), VariantAnnotation::alt(sample), "]",
                                    as.character(Biostrings::subseq(sample$context, start = (context-1)/2 + 2, end = context)))

  # Remove context to save memory
  sample$context <- NULL

  return(sample)
}
