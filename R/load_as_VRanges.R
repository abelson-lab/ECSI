#' Load as VRanges
#'
#' Load varscan pileup2cns output from a given sample as a \code{VRanges} object
#'
#' @param sample_name Name of the sample
#' @param sample_path Sample file location
#' @param genome Reference genome to use, default is hg19
# @param MAPQ_cutoff Minimum acceptable MAPQ score for variants to keep
#' @param metadata Logical. Whether to include metadata (VAF, quality scores, strand-specific counts)
# @importFrom data.table fread
# @importFrom methods as
#' @importFrom dplyr "%>%"
# @importMethodsFrom S4Vectors mcols
# @importMethodsFrom GenomeInfoDb genome
# @importClassesFrom VariantAnnotation VRanges
#' @export
#' @examples
#' \dontrun{
#' variants <- load_as_VRanges(sample_name = "pt123", sample_path = "./patient_123_pileup2cns",
#' genome = "hg19", MAPQ_cutoff = 59, metadata = TRUE)
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

### TO DO LIST:
#	- Check input to ensure it is suitable
#	- Accommodate hg38 as well

load_as_VRanges <-
function(sample_name, sample_path, genome = "hg19", metadata = TRUE) {

  # Load in as dataframe
  varscan_output_df <- data.table::fread(sample_path) %>%
    dplyr::filter(Reads2 != 0)

  # Convert to VRanges
  varscan_output <- with(varscan_output_df, VariantAnnotation::VRanges(
    seqnames = paste0("chr",Chrom),
    ranges = IRanges(Position, Position),
    ref = Ref, alt = VarAllele,
    refDepth = Reads1, altDepth = Reads2, sampleNames = sample_name))

  if(metadata==TRUE) {

    # Add metadata
    S4Vectors::mcols(varscan_output) <- varscan_output_df %>%
      dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
      dplyr::select(VAF, Qual1, Qual2, MapQual1, MapQual2, Reads1Plus, Reads1Minus, Reads2Plus, Reads2Minus)

    # Convert quality scores to Rle to save memory
    varscan_output$Qual1 <- methods::as(varscan_output$Qual1, "Rle")
    varscan_output$Qual2 <- methods::as(varscan_output$Qual2, "Rle")
    varscan_output$MapQual1 <- methods::as(varscan_output$MapQual1, "Rle")
    varscan_output$MapQual2 <- methods::as(varscan_output$MapQual2, "Rle")

  } else if (metadata == "VAF"){

    # Add metadata
    S4Vectors::mcols(varscan_output) <- varscan_output_df %>%
      dplyr::mutate(VAF = Reads2 / (Reads1 + Reads2)) %>%
      dplyr::select(VAF)
  }

  # specify genome
  GenomeInfoDb::genome(varscan_output) = genome    # for now, must be hg19

  # clean up and remove dataframe
  rm(varscan_output_df)

  return(varscan_output)
}
