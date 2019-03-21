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
#' @param recurrent_mutations \code{VRanges} object with chr, pos ref, alt of frequently mutated alleles to remove from model input. \code{GRanges} objects are also accepted, in which case filtering will occur by position.
#' @export
#' @examples
#' \dontrun{
#' # Get flagged alleles and cosmic mutations
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "./heme_COSMIC.csv")
#' flagged_alleles <- get_flagged_alleles(all_sample_names, all_sample_paths, exclude_cosmic_mutations = TRUE,
#'     cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
#'
#' # Load and annotate sample
#' samp <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' samp <- sequence_context(samp)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' annotated_samp <- annotate_MAF(varscan_output = variants,
#'     MAF_database = MafDb.gnomADex.r2.1.hs37d5, genome = "hg19")
#'
#' # Filter model input
#' samp_model_input <- filter_model_input(model_input = annotated_samp, flagged_alleles = flagged_alleles,
#'     filter_cosmic_mutations = TRUE, cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 10)
#' }
#' @return This function returns a filtered \code{VRanges} object.
#'
#'

filter_model_input <-
function(model_input, flagged_alleles = NA, MAF_cutoff = 0.001, VAF_cutoff = 0.05, filter_cosmic_mutations = FALSE, cosmic_mutations, cosmic_mut_frequency = 10,
         MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59, recurrent_mutations = NA) {

  # Check that input is VRanges
  if(class(model_input) != "VRanges"){ stop('Input "model_input" needs to be VRanges object') }

  # keep variants below MAF cutoff (keep NA because that means 0)
  model_input <- model_input[which(model_input$AF < MAF_cutoff | is.na(model_input$AF)),]
  # keep variants below VAF cutoff
  model_input <- model_input[which(model_input$VAF < VAF_cutoff),]

  # filter MAPQ
  model_input <- filter_MAPQ(model_input, MAPQ_cutoff_ref = MAPQ_cutoff_ref, MAPQ_cutoff_alt = MAPQ_cutoff_alt)

  # If user provides recurrent mutations
  if(class(flagged_alleles) %in% c("VRanges", "GRanges")){
    # remove variants that overlap with cosmic mutations
    model_input <- subtract_VRanges(model_input, flagged_alleles)
  }

  # if filtering cosmic mutations
  if(filter_cosmic_mutations == TRUE){
    # filter for frequency above 10
    cosmic_mutations <- cosmic_mutations[cosmic_mutations$hemCOSMIC_DC >= cosmic_mut_frequency]
    # remove variants that overlap with cosmic mutations
    model_input <- subtract_VRanges(model_input, cosmic_mutations)
  }

  # If user provides recurrent mutations
  if(class(recurrent_mutations) %in% c("VRanges", "GRanges")){
    # remove variants that overlap with cosmic mutations
    model_input <- subtract_VRanges(model_input, recurrent_mutations)
  }

  return(model_input)
}
