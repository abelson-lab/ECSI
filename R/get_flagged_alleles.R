#' Get Flagged Alleles
#'
#' Flag alleles that are present in too many samples at high variant allele frequencies as potential errors.
#'
#' @param files List of paths to all variant files from varscan pileup2cns output
#' @param exclude_cosmic_mutations Logical indicating whether or not to exclude cosmic mutations from flagged SNPs
#' @param cosmic_mutations \code{VRanges} object with position and substitution of excluded cosmic mutations
#' @param cosmic_mut_frequency Mutations with this frequency or above in the cosmic database will be excluded
# @importClassesFrom VariantAnnotation VRanges SimpleVRangesList
# @importMethodsFrom S4Vectors mcols queryHits
# @importMethodsFrom GenomicRanges findOverlaps
# @importMethodsFrom BiocGenerics match
#' @export
#' @examples
#' \dontrun{
#' files <- list.files(path = "./", pattern = "sample")
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "./heme_COSMIC.csv")
#'
#' flagged_alleles <- get_flagged_alleles(files, exclude_cosmic_mutations = TRUE,
#' cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
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
#	- Apply basic quality filter to input prior to flagging variants
# 	- Parameter to tune quantile starting point
#	- Parameter to tune fisher test cut-off

# may have to input functions as arguments

get_flagged_alleles <-
function(files, exclude_cosmic_mutations = FALSE, cosmic_mutations, cosmic_mut_frequency = 3){

  ### GRANGES LIST TO STORE ALL THE SAMPLES BY POSITION AND VAF
  vrlist <- VariantAnnotation::VRangesList()

  for(samp_path in files){
    # samp name
    samp_name <- substr(samp_path,4,12)
    # get sample as VRanges and annotate with sequence context and MAF
    samp <- load_as_VRanges(samp_name, samp_path, metadata = "VAF")
    # save as item in grlist
    vrlist[[samp_name]] <- samp
  }

  ### CREATE THE VARIANTS MATRIX TO FEED INTO TAGGED SNPS

  # all unique positions
  variants = unique(unlist(vrlist))
  S4Vectors::mcols(variants)$VAF <- NULL

  for (i in 1:length(vrlist)){
    samp_name <- names(vrlist)[i]
    samp <- vrlist[[i]]
    S4Vectors::mcols(variants)[BiocGenerics::match(samp, variants), samp_name] = samp$VAF
  }

  ### TAG THE FREQUENT SNPS
  flagged_alleles <- flag_alleles(variants)

  # exclude the cosmic mutations
  if(exclude_cosmic_mutations == TRUE){
    # filter for frequency above 3
    cosmic_mutations <- cosmic_mutations[cosmic_mutations$hemCOSMIC_DC >= cosmic_mut_frequency]
    # removed flagged alleles overlapping with cosmic mutations
    flagged_alleles <- setdiff_VRanges(flagged_alleles, cosmic_mutations)
  }

  return(flagged_alleles)
}
