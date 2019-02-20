#' Annotate MAF
#'
#' Annotate VRanges object with population minor allele frequency from gnomAD
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param MAF_database \code{GScores} object with MAFs for each position
#' @param genome Genome of input object, accepts either hg19 or hg38
#' @export
#' @examples
#' \dontrun{
#' # Sample with hg19 reference, annotated with gnomAD exome hs37d5
#' variants <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' variants <- annotate_MAF(varscan_output = variants,
#'     MAF_database = MafDb.gnomADex.r2.1.hs37d5, genome = "hg19")
#' }
#' \dontrun{
#' # Sample with hg38 reference, annotated with gnomAD genome GRCh38
#' variants <- load_as_VRanges(sample_name = "pt321",
#'     sample_path = "./patient_321_pileup2cns", genome = "hg38", metadata = TRUE)
#' library(MafDb.gnomAD.r2.0.1.GRCh38)
#' variants <- annotate_MAF(varscan_output = variants,
#'     MAF_database = MafDb.gnomAD.r2.0.1.GRCh38, genome = "hg38")
#' }
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{MAF}.

### TO DO LIST:
# Check that varscan_output is GRanges or VRanges Object
# Check that MAF_database is a GScores object

annotate_MAF <-
  function(varscan_output, MAF_database, genome = c("hg19", "hg38")){

    # temp convert genome name to match gscores
    if(genome == "hg19"){
      GenomeInfoDb::genome(varscan_output) <- "hs37d5"
    } else if(genome == "hg38"){
      GenomeInfoDb::genome(varscan_output) <- "GRCh38"
    }

    # Annotate MAF
    varscan_output_anno <- GenomicScores::gscores(MAF_database, varscan_output)
    GenomeInfoDb::seqlevelsStyle(varscan_output_anno) <- "UCSC"

    # convert genome name back to input
    if(genome == "hg19"){
      GenomeInfoDb::genome(varscan_output_anno) <- "hg19"
    } else if(genome == "hg38"){
      GenomeInfoDb::genome(varscan_output_anno) <- "hg38"
    }

    return(varscan_output_anno)
  }
