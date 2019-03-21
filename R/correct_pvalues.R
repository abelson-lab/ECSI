#' p-value correction
#'
#' Apply bonferroni multiple testing correction to model p-values for each position.
#' For mutation calling in a narrow genomic range of specified hotspots, we recommend correcting across all samples.
#' For mutation calling across all genomic positions sequenced, we recommend applying correction by sample.
#'
#' @param variant_calls \code{VRanges} object from call_all_variants output
#' @param by_sample Logical. Option to apply multiple testing correction by sample.
#' @export
#' @examples
#' \dontrun{
#' variants_corrected <- correct_pvalue(variant_calls = variants, by_sample = FALSE)
#' }
#' @return This function returns a \code{VRanges} object with a new metadata column: corrected_pvalue

correct_pvalues <-
function(variant_calls, by_sample = FALSE){

  if(by_sample == TRUE){
    # get number of p values generated for each sample
    sample_counts <- table(sampleNames(variant_calls))

    # get corrected pvalues
    variant_calls$corrected_pvalue <- variant_calls$model_Pvalue * as.numeric(sample_counts[as.character(sampleNames(variant_calls))])
    variant_calls$corrected_pvalue <- ifelse(variant_calls$corrected_pvalue > 1, 1, variant_calls$corrected_pvalue)

  } else if(by_sample == FALSE){
    # get corrected pvalues
    variant_calls$corrected_pvalue <- variant_calls$model_Pvalue * length(variant_calls)
    variant_calls$corrected_pvalue <- ifelse(variant_calls$corrected_pvalue > 1, 1, variant_calls$corrected_pvalue)

  }

  return(variant_calls)
}
