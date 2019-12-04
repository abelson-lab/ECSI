#' Load File with Recurrent Mutations to exclude
#'
#' Loads in recurrent mutations from tab or comma-delimited file with four columns (Chr, Pos, Ref, Alt) and no header, presenting it as a VRanges file.
#'
#' @param file Location of file with recurrent mutations
#' @param genome Reference genome to use, default is hg19
#' @export
#' @examples
#' \dontrun{
#' recurrent <- load_recurrent_mutations(file = "COSMIC_Heme_mutations_freq10.txt", genome = "hg19")
#' }
#' @return This function returns a \code{VRanges} object including:
#' \itemize{
#'	\item seqnames
#' 	\item ranges
#' 	\item ref
#' 	\item alt
#' }

load_recurrent_mutations <-
function(file, genome = c("hg19", "hg38")) {

  # Check that genome was specified by the user
  if(length(genome) > 1){ stop("Need to specify genome, either hg19 or hg38.") }

  # load as datatable (automatic handling of headers and tab/comma delimiter)
  recurrent <- data.table::fread(file)

  # convert to data.frame, fill in colnames, standardize chr names
  recurrent_df <- as.data.frame(recurrent)
  colnames(recurrent_df)[1:4] <- c("chr", "pos", "ref", "alt")
  recurrent_df$chr <- ifelse(stringr::str_detect(recurrent_df$chr, "[Cc][Hh][Rr]"), stringr::str_replace(recurrent_df$chr, "[Cc][Hh][Rr]", ""), recurrent_df$chr)

  # Convert to VRanges
  recurrent_vr <- VariantAnnotation::VRanges(
    seqnames = paste0("chr",recurrent_df$chr),
    ranges = IRanges::IRanges(recurrent_df$pos, recurrent_df$pos),
    ref = recurrent_df$ref,
    alt = recurrent_df$alt)

  # specify genome
  GenomeInfoDb::genome(recurrent_vr) = genome

  return(recurrent_vr)
}
