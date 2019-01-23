#' Annotate MAF
#'
#' Annotate VRanges object with population minor allele frequency from gnomAD
#'
#' @param varscan_output \code{VRanges} object of the varscan pileup2cns output
#' @param MAF_database \code{GScores} object with MAFs for each position
# @importClassesFrom VariantAnnotation VRanges
# @importFrom MafDb.gnomAD.r2.0.1.GRCh38 MafDb.gnomAD.r2.0.1.GRCh38
# @importFrom GenomicScores gscores
# @importFrom utils data
# @importFrom GenomeInfoDb genome
# @importFrom rtracklayer liftOver
#' @export
#' @examples
#' \dontrun{
#' variants <- load_as_VRanges(sample_name = "pt123",
#' sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' variants <- annotate_MAF(varscan_output = variants, MAF_database = MafDb.gnomADex.r2.1.hs37d5)
#' }
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
  function(varscan_output, MAF_database){

#     # Liftover from hg19 to hg38
#     #hg19tohg38 <- import.chain(liftOver_chain_hg19toHg38)
#     utils::data("hg19tohg38", envir = environment())
#     varscan_output_hg38 <- unlist(rtracklayer::liftOver(varscan_output, hg19tohg38))
#     GenomeInfoDb::genome(varscan_output_hg38) <- "GRCh38"
#
#     # Annotate with MAF from Exac
#     varscan_output_hg38_anno <- GenomicScores::gscores(MafDb.gnomAD.r2.0.1.GRCh38::MafDb.gnomAD.r2.0.1.GRCh38, varscan_output_hg38)
#     GenomeInfoDb::seqlevelsStyle(varscan_output_hg38_anno) <- "UCSC"
#
#     # Liftover from hg38 back to hg19
#     #hg38tohg19 <- import.chain(liftOver_chain_hg38toHg19)
#     utils::data("hg38tohg19", envir = environment())
#     varscan_output_anno <- unlist(rtracklayer::liftOver(varscan_output_hg38_anno, hg38tohg19))
#     GenomeInfoDb::genome(varscan_output_anno) <- "hg19"
#
#     rm(hg19tohg38, hg38tohg19)

    # Annotate MAF
    varscan_output_anno <- GenomicScores::gscores(MAF_database, varscan_output)
    GenomeInfoDb::seqlevelsStyle(varscan_output_anno) <- "UCSC"

    return(varscan_output_anno)
  }
