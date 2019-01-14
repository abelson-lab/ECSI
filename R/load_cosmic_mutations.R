#' Load Cosmic Mutations
#'
#' Loads in Cosmic mutations from file and present it as a VRanges object.
#'
#' @param cosmic_mutations_path Path to csv file containing COSMIC mutations
# @importFrom dplyr select
# @importFrom data.table fread
# @importMethodsFrom S4Vectors mcols
# @importClassesFrom IRanges IRanges
# @importClassesFrom VariantAnnotation VRanges
#' @export
#' @examples
#' \dontrun{
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "heme_COSMIC.csv")
#' }
#' @return This function returns a \code{VRanges} object including:
#' \itemize{
#'	\item seqnames
#' 	\item ranges
#'	\item ref
#'	\item alt
#'	\item metadata (including observed frequency)
#' }

load_cosmic_mutations <-
function(cosmic_mutations_path = "hemCOSMIC.csv") {
  # load in as dataframe
  cosmic_mutations_df <- data.table::fread(cosmic_mutations_path)

  # convert to VRanges
  cosmic_mutations <- with(cosmic_mutations_df, VariantAnnotation::VRanges(
    seqnames = paste0("chr",Chr), ranges = IRanges::IRanges(Start, End), ref = Ref, alt = Alt))

  # load in metadata
  S4Vectors::mcols(cosmic_mutations) <- dplyr::select(cosmic_mutations_df, -c("Chr", "Start", "End", "Ref", "Alt"))

  return(cosmic_mutations)
}
