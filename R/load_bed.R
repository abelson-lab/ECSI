#' Load Bed File as GRanges
#'
#' Loads in sequence ranges from bed file and presents it as a GRanges object. Can handle headers and different Chr naming conventions.
#'
#' @param bed_path Path to bed file
#' @param genome Reference genome to use, default is hg19
#' @export
#' @examples
#' \dontrun{
#' regions_of_interest <- load_bed(bed_path = "AML_exons.bed", genome = "hg19")
#' }
#' @return This function returns a \code{GRanges} object including:
#' \itemize{
#'	\item seqnames
#' 	\item ranges
#' }

load_bed <-
function(bed_path, genome = "hg19") {
  # load as datatable (automatic handling of headers)
  bed_dt <- data.table::fread(bedpath)

  # convert to data.frame, fill in colnames, standardize chr names
  bed_df <- as.data.frame(bed_dt)
  colnames(bed_df)[1:3] <- c("chr", "start", "stop")
  bed_df <- dplyr::mutate(bed_df, chr = ifelse(stringr::str_detect(chr, "[Cc][Hh][Rr]"), stringr::str_replace(chr, "[Cc][Hh][Rr]", ""), chr))

  # convert to GRanges
  bed_gr <- GenomicRanges::GRanges(
    seqnames = paste0("chr", bed_df$chr),
    ranges = IRanges::IRanges(bed_df$start, bed_df$stop))

  # specify genome
  GenomeInfoDb::genome(bed_gr) = genome

  return(bed_gr)
}
