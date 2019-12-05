#' Filter Model Input
#'
#' Filter variants to remove flagged alleles, polymorphisms, cosmic mutations, and high VAF prior to error model generation.
#'
#' @param model_input \code{VRanges} object annotated with mutation context and population minor allele frequency
#' @param flagged_alleles \code{VRanges} object with high VAF alleles flagged as being present in too many samples
#' @param MAF_cutoff Population Minor Allele Frequency cutoff: variants at or above this cutoff are excluded. Default is 0.001
#' @param VAF_cutoff Sample Variant Allele Frequency cutoff: variants at or above this cutoff are excluded. Default is 0.05
#' @param MAPQ_cutoff Minimum acceptable MAPQ score; positions below this cutoff will be excluded. Default is 59
#' @param recurrent_mutations \code{VRanges} object with chr, pos ref, alt of frequently mutated alleles to remove from model input. \code{GRanges} objects are also accepted, in which case filtering will occur by position.
#' @export
#' @examples
#' \dontrun{
#' # Get flagged alleles and cosmic mutations
#' hemeCOSMIC_10 <- load_recurrent_mutations("example_data/COSMIC_heme_freq10.txt", genome = "hg19")
#' flagged_alleles <- get_flagged_alleles(all_sample_names, all_sample_paths,
#'     exclude_cosmic_mutations = TRUE, cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
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
#' samp_model_input <- filter_model_input(model_input = annotated_samp,
#'     flagged_alleles = flagged_alleles, recurrent_mutations = hemeCOSMIC_10)
#' }
#' @return This function returns a filtered \code{VRanges} object.


filter_model_input <-
function(model_input, flagged_alleles = NA, MAF_cutoff = 0.001, VAF_cutoff = 0.05, MAPQ_cutoff = 59, recurrent_mutations = NA) {

  # Check that input is VRanges
  if(class(model_input) != "VRanges"){ stop('Input "model_input" needs to be VRanges object') }

  # keep variants below MAF cutoff (keep NA because that means 0)
  model_input <- model_input[which(model_input$MAF < MAF_cutoff | is.na(model_input$MAF)),]
  # keep variants below VAF cutoff
  model_input <- model_input[which(model_input$VAF < VAF_cutoff),]

  # filter MAPQ
  model_input <- filter_MAPQ(model_input, MAPQ_cutoff = MAPQ_cutoff)

  # If user provides flagged alleles
  if(class(flagged_alleles) %in% c("VRanges", "GRanges")){
    # remove variants that overlap with cosmic mutations
    model_input <- subtract_VRanges(model_input, flagged_alleles)
  }

  # If user provides recurrent mutations
  if(class(recurrent_mutations) %in% c("VRanges", "GRanges")){
    # remove variants that overlap with cosmic mutations
    model_input <- subtract_VRanges(model_input, recurrent_mutations)
  }

  return(model_input)
}
