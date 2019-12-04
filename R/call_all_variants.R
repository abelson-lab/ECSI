#' Call All Variants
#'
#' Call all variants from the varscan pileup2cns output and the context-specific error models
#' In the rare case that no model exists for that context (not enough alternate alleles), then this function uses the varscan assigned pvalue
#'
#' @param sample \code{VRanges} object of the varscan pileup2cns output annotated with variant context
#' @param samp_models \code{Data.Frame} of context specific error models generated from \code{generate_all_models}
#' @importFrom foreach "%dopar%"
#' @export
#' @examples
#' \dontrun{
#' # Get flagged alleles and cosmic mutations
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "./heme_COSMIC.csv")
#' flagged_alleles <- get_flagged_alleles(all_sample_names, all_sample_paths,
#'     exclude_cosmic_mutations = TRUE,cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
#'
#' # Load and annotate sample
#' sample <- load_as_VRanges(sample_name = "pt123",
#'     sample_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' sample <- sequence_context(sample)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' annotated_samp <- annotate_MAF(varscan_output = variants,
#'     MAF_database = MafDb.gnomADex.r2.1.hs37d5, genome = "hg19")
#'
#' # Filter model input
#' samp_model_input <- filter_model_input(model_input = annotated_samp,
#'     flagged_alleles = flagged_alleles, filter_cosmic_mutations = TRUE,
#'     cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 10)
#'
#' # Generate the error models for this sample
#' samp_models <- generate_all_models(sample = samp_model_input, plots = FALSE)
#'
#' # Call variants for the sample based on the previously generated error models
#' variant_calls <- call_all_variants(annotated_samp, samp_models)
#' }
#' @return This function returns a \code{VRanges} with the following metadata:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used
#'	\item Model pvalue
#'	}

call_all_variants <-
function(sample, samp_models) {

  # Check that input is VRanges
  if(class(sample) != "VRanges"){ stop('Input "sample" needs to be VRanges object') }

  # Create model and model_Pvalue columns
  sample$model <- NA
  sample$model_Pvalue <- NA

  # Get 192 unique trinucleotide combos
  flanking_seqs <- unique(sample$FlankingSeqGroup)

  # for each signature in trinucleotide combo
    # compare alternate allele count with
  m <- foreach::foreach(i=1:length(flanking_seqs), .combine='c') %dopar% {
    call_variants(i, sample, flanking_seqs, samp_models)
  }

  # index for all the variant calls
  m1=unlist(m[c(TRUE, FALSE, FALSE)])
  # pvals corresponding to variant call in each entry of m1
  m2=unlist(m[c(FALSE, TRUE, FALSE)])
  # model used corresponding to variant call in each entry of m1
  m3=unlist(m[c(FALSE, FALSE, TRUE)])
  # assign model used in original dataframe
  sample$model[m1]=m3
  # assign pvalue in original dataframe
  sample$model_Pvalue[m1]=m2
  # NOTE THAT THIS IS BEFORE MULTIPLE TESTING

  # if no error model available for the context, use the varscan pvalue
  sample[sample$model == "None"]$model_Pvalue = sample[sample$model == "None"]$Varscan_Pval
  sample$Varscan_Pval <- NULL

  # covert sample$model to RLE to save memory
  sample$model <- methods::as(sample$model, "Rle")

  return(sample)
}
