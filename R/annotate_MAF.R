#' Annotate MAF
#'
#' Annotate VRanges object with population minor allele frequency from gnomAD
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @import VariantAnnotation
#' @import MafDb.gnomAD.r2.0.1.GRCh38
#' @import GenomicScores
#' @export
#' @examples
#' variants <- load_as_VRanges(sample_name = "pt123", sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' variants <- annotate_MAF(varscan_output = variants)
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{MAF}.

### TO DO LIST:
# Check that input is GRanges Object
# Check that input is hg38 and adapt accordingly
# Load liftover chains as part of package
# close connections warning with chain files, need to resolve this

# can probably resolve by loading chain files directly into package workspace and instead of loading within the function


### Assuming chain files are already loaded in the env
# function(varscan_output, liftOver_chain_hg19toHg38 = "hg19ToHg38.over.chain", liftOver_chain_hg38toHg19 = "hg38ToHg19.over.chain"){

annotate_MAF <-
function(varscan_output){

  # Liftover from hg19 to hg38
  #hg19tohg38 <- import.chain(liftOver_chain_hg19toHg38)
  varscan_output_hg38 <- unlist(liftOver(varscan_output, hg19tohg38))
  genome(varscan_output_hg38) <- "GRCh38"

  # Annotate with MAF from Exac
  varscan_output_hg38_anno <- gscores(MafDb.gnomAD.r2.0.1.GRCh38, varscan_output_hg38)
  seqlevelsStyle(varscan_output_hg38_anno) <- "UCSC"

  # Liftover from hg38 back to hg19
  #hg38tohg19 <- import.chain(liftOver_chain_hg38toHg19)
  varscan_output_anno <- unlist(liftOver(varscan_output_hg38_anno, hg38tohg19))
  genome(varscan_output_anno) <- "hg19"

  return(varscan_output_anno)
}
