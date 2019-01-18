#' Call All Variants
#'
#' Call all variants from the varscan pileup2cns output and the context-specific error models
#'
#' @param samp \code{VRanges} object of the varscan pileup2cns output annotated with variant context
#' @param samp_models \code{Data.Frame} of context specific error models generated from \code{generate_all_models}
# @importClassesFrom VariantAnnotation VRanges
#' @importFrom foreach "%dopar%"
# @import doParallel
#' @export
#' @examples
#' \dontrun{
#' model_input <- filter_model_input(model_input = samp,
#' flagged_alleles = flagged_alleles, filter_cosmic_mutations = TRUE,
#' cosmic_mutations = hemeCOSMIC, cosmic_mut_frequency = 10)
#' samp_models <- generate_all_models(samp = model_input, plots = FALSE)
#' variant_calls <- call_all_variants(samp, samp_models)
#' }
#' @return This function returns a \code{VRanges} with the following metadata:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used
#'	\item Model pvalue
#'	}

### TO DO LIST:
#	- Check input to ensure it is suitable
# - Is data.table even needed

# may have to input functions as arguments

call_all_variants <-
function(samp, samp_models) {

  # Create model and model_Pvalue columns
  samp$model <- NA
  samp$model_Pvalue <- NA

  # Get 192 unique trinucleotide combos
  flanking_seqs <- unique(samp$FlankingSeqGroup)

  # for each signature in trinucleotide combo
    # compare alternate allele count with
  m <- foreach::foreach(i=1:length(flanking_seqs), .combine='c') %dopar% {
    call_variants(i, samp, flanking_seqs, samp_models)
  }

  # index for all the variant calls
  m1=unlist(m[c(TRUE, FALSE, FALSE)])
  # pvals corresponding to variant call in each entry of m1
  m2=unlist(m[c(FALSE, TRUE, FALSE)])
  # model used corresponding to variant call in each entry of m1
  m3=unlist(m[c(FALSE, FALSE, TRUE)])
  # assign model used in original dataframe
  samp$model[m1]=m3
  # assign pvalue in original dataframe
  samp$model_Pvalue[m1]=m2
    # NOTE THAT THIS IS BEFORE MULTIPLE TESTING
    # NOTE THAT IF MODEL IS NOT AVAILABLE, PVALUE IS 1. THIS IS NOT IDEAL.

  # covert samp$model to RLE to save memory
  samp$model <- methods::as(samp$model, "Rle")

  return(samp)
}
