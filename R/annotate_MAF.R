#' Annotate MAF
#'
#' Annotate VRanges object with population minor allele frequency from gnomAD
#'
#' @param sample \code{VRanges} object of the varscan pileup2cns output
#' @param MAF_database \code{GScores} object with MAFs for each position
#' @param genome Genome of input object, accepts either hg19 or hg38
#' @export
#' @examples
#' \dontrun{
#' # Sample with hg19 reference, annotated with gnomAD exome hs37d5
#' variants <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' variants <- annotate_MAF(sample = variants,
#'     MAF_database = MafDb.gnomADex.r2.1.hs37d5, genome = "hg19")
#' }
#' \dontrun{
#' # Sample with hg38 reference, annotated with gnomAD genome GRCh38
#' variants <- load_as_VRanges(sample_name = "pt321",
#'     sample_path = "./patient_321_pileup2cns", genome = "hg38", metadata = TRUE)
#' library(MafDb.gnomAD.r2.0.1.GRCh38)
#' variants <- annotate_MAF(sample = variants,
#'     MAF_database = MafDb.gnomAD.r2.0.1.GRCh38, genome = "hg38")
#' }
#' @return This function returns the \code{VRanges} object with the additional metadata column \code{MAF}.

### TO DO LIST:
# Check that sample is GRanges or VRanges Object
# Check that MAF_database is a GScores object

annotate_MAF <-
  function(sample, MAF_database, genome = c("hg19", "hg38")){

    # Check that genome was specified by the user
    if(length(genome) > 1){ stop("Need to specify genome, either hg19 or hg38.") }
    # Check that input is VRanges
    if(class(sample) != "VRanges"){ stop('Input "sample" needs to be VRanges object') }

    # temp convert genome name to match gscores
    if(genome == "hg19"){
      GenomeInfoDb::genome(sample) <- "hs37d5"
    } else if(genome == "hg38"){
      GenomeInfoDb::genome(sample) <- "GRCh38"
    } else {
      stop(paste("Genome",genome,"not recognized. Genome must be either hg19 or hg38."))
    }

    # Annotate MAF
    sample_anno <- GenomicScores::gscores(MAF_database, sample)
    GenomeInfoDb::seqlevelsStyle(sample_anno) <- "UCSC"
    sample_anno$MAF <- sample_anno$AF; sample_anno$AF <- NULL

    # convert genome name back to input
    if(genome == "hg19"){
      GenomeInfoDb::genome(sample_anno) <- "hg19"
    } else if(genome == "hg38"){
      GenomeInfoDb::genome(sample_anno) <- "hg38"
    }

    return(sample_anno)
  }
