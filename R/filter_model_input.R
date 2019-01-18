#' Filter Model Input
#'
#' Filter variants to remove flagged alleles, polymorphisms, cosmic mutations, and high VAF prior to error model generation.
#'
#' @param model_input \code{VRanges} object annotated with mutation context and population minor allele frequency
#' @param flagged_alleles \code{VRanges} object with high VAF alleles flagged as being present in too many samples
#' @param MAF_cutoff Population Minor Allele Frequency cutoff: variants at or above this cutoff are excluded
#' @param VAF_cutoff Sample Variant Allele Frequency cutoff: variants at or above this cutoff are excluded
#' @param filter_cosmic_mutations Logical indicating whether or not to filter cosmic mutations
#' @param cosmic_mutations \code{VRanges} object with position and substitution of excluded cosmic mutations
#' @param cosmic_mut_frequency Mutations with this frequency or above in the cosmic database will be excluded
#' @param MAPQ_cutoff_ref Minimum acceptable MAPQ score for reference allele
#' @param MAPQ_cutoff_alt Minimum acceptable MAPQ score for alternate allele
#' @param filter_custom_alleles Logical. Whether to filter out custom alleles (e.g. specific mutational sites) from model input
#' @param custom_allele_positions \code{VRanges} object with chr, pos, ref, alt of custom alleles to filter out of model input
# @importMethodsFrom S4Vectors queryHits
# @importMethodsFrom GenomicRanges findOverlaps
# @importClassesFrom VariantAnnotation VRanges
#' @export
#' @examples
#' \dontrun{
#' samp_models <- filter_model_input(model_input = annotated_samp, flagged_alleles = flagged_alleles,
#'  filter_cosmic_mutations = TRUE, cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 10)
#' }
#' @return This function returns a filtered \code{VRanges} object.

### TO DO LIST:
#	- Check input to ensure it is suitable
#	- Provide option for flagged alleles

# may have to input functions as arguments

filter_model_input <-
function(model_input, flagged_alleles, MAF_cutoff = 0.001, VAF_cutoff = 0.05, filter_cosmic_mutations = FALSE, cosmic_mutations, cosmic_mut_frequency = 10,
         MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59, filter_custom_alleles = FALSE, custom_allele_positions) {

  # keep variants below MAF cutoff
  model_input <- model_input[which(model_input$AF < MAF_cutoff),]
  # keep variants below VAF cutoff
  model_input <- model_input[which(model_input$VAF < VAF_cutoff),]

  # filter MAPQ
  model_input <- filter_MAPQ[model_input, MAPQ_cutoff_ref = MAPQ_cutoff_ref, MAPQ_cutoff_alt = MAPQ_cutoff_alt]

  # remove variants that overlap with flagged alleles
  model_input <- model_input[-S4Vectors::queryHits(GenomicRanges::findOverlaps(model_input, flagged_alleles, type = "any"))

  # if filtering cosmic mutations
  if(filter_cosmic_mutations == TRUE){
    # filter for frequency above 10
    cosmic_mutations <- cosmic_mutations[cosmic_mutations$hemCOSMIC_DC >= cosmic_mut_frequency]
    # remove variants that overlap with cosmic mutations
    model_input <- model_input[-S4Vectors::queryHits(GenomicRanges::findOverlaps(model_input, cosmic_mutations, type = "any"))]
  }

  # if removing custom alleles
  if(filter_custom_alleles == TRUE){
    # remove variants that overlap with cosmic mutations
    model_input <- model_input[-S4Vectors::queryHits(GenomicRanges::findOverlaps(model_input, custom_allele_positions, type = "any"))]
  }

  return(model_input)
}
