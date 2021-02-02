#' Load as VRanges
#'
#' Load varscan pileup2cns output from a given sample as a \code{VRanges} object
#'
#' @param sample_name Name of the sample
#' @param sample_path Sample file location
#' @param genome Reference genome to use.
#' @param metadata Logical. Whether to include metadata (VAF, quality scores, strand-specific counts)
#' @importFrom dplyr "%>%"
#' @export
#' @examples
#' \dontrun{
#' variants <- load_as_VRanges(sample_name = "pt123", sample_path = "./patient_123_pileup2cns",
#' genome = "hg19", metadata = TRUE)
#' }
#' @return This function returns a \code{VRanges} object with the following information:
#' \itemize{
#'	\item seqnames
#' 	\item ranges
#'	\item ref
#'	\item alt
#'	\item refdepth
#'	\item altdepth
#'	\item sampleNames
#'	\item metadata (optional)
#'	}

load_as_VRanges <-
function(sample_name, sample_path, genome, metadata = TRUE) {

  # Load in as dataframe
  varscan_df <- data.table::fread(sample_path) %>%
    dplyr::filter(Reads2 != 0)

  if(metadata == "VAF"){

  # Convert to VRanges but without refDepth or altDepth
  varscan_output <- VariantAnnotation::VRanges(
    if(length(grep("chr",varscan_df$Chrom))==0) {seqnames = paste0("chr",varscan_df$Chrom)} else {seqnames = varscan_df$Chrom} ,
    ranges = IRanges::IRanges(varscan_df$Position, varscan_df$Position),
    ref = varscan_df$Ref,
    alt = varscan_df$VarAllele,
    sampleNames = sample_name)

  # Add metadata
  S4Vectors::mcols(varscan_output) <- varscan_df %>%
    dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
    dplyr::select(VAF)

  } else {

    # Convert to VRanges
    varscan_output <- VariantAnnotation::VRanges(
      if(length(grep("chr",varscan_df$Chrom))==0) {seqnames = paste0("chr",varscan_df$Chrom)} else {seqnames = varscan_df$Chrom} ,
      ranges = IRanges::IRanges(varscan_df$Position, varscan_df$Position),
      ref = varscan_df$Ref,
      alt = varscan_df$VarAllele,
      refDepth = varscan_df$Reads1,
      altDepth = varscan_df$Reads2,
      sampleNames = sample_name)

    if(metadata==TRUE) {

      # Add metadata
      S4Vectors::mcols(varscan_output) <- varscan_df %>%
        dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
        dplyr::select(VAF, Qual1, Qual2, MapQual1, MapQual2, Reads1Plus, Reads1Minus, Reads2Plus, Reads2Minus, "Varscan_Pval" = Pvalue)

      # Convert quality scores to Rle to save memory
      varscan_output$Qual1 <- methods::as(varscan_output$Qual1, "Rle")
      varscan_output$Qual2 <- methods::as(varscan_output$Qual2, "Rle")
      varscan_output$MapQual1 <- methods::as(varscan_output$MapQual1, "Rle")
      varscan_output$MapQual2 <- methods::as(varscan_output$MapQual2, "Rle")

    }
  }

  # specify genome
  GenomeInfoDb::genome(varscan_output) = genome

  # clean up and remove dataframe
  rm(varscan_df)

  return(varscan_output)
}
