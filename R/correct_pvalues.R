#' p-value correction
#'
#' Apply multiple testing correction (FDR or Bonferroni) to model p-values for each position.
#' If input is VRangesList, then pvalue correction will be done on a sample-by-sample basis. If input is VRanges, then pvalue correction will be done across all samples.
#' For mutation calling in a narrow genomic range of specified hotspots, we recommend correcting across all samples.
#' For mutation calling across all genomic positions sequenced, we recommend applying correction by sample.
#'
#' @param variant_calls \code{VRangesList} or \code{VRanges} object from call_all_variants output
#' @param method pvalue correction method to pass onto p.adjust
#' @export
#' @examples
#' \dontrun{
#' variants_corrected <- correct_pvalue(variant_calls = variants, method = 'fdr')
#' }
#' @return This function returns a \code{VRangesList} or \code{VRanges} object with a new metadata column: corrected_pvalue

correct_pvalues <-
function(variant_calls, method = 'fdr'){

  # if VRangesList, correct pvalues by sample
  if(class(variant_calls) %in% c("VRangesList", "SimpleVRangesList")){
    variant_calls <- S4Vectors::endoapply(variant_calls, function(x){x$corrected_pvalue <- p.adjust(x$model_Pvalue, method); return(x)})

  } else if(class(variant_calls) %in% c("VRanges")){
    variant_calls$corrected_pvalue <- p.adjust(x$model_Pvalue, method)

  } else {
    stop('input "variant_calls" needs to be of class "VRanges" or "VRangesList"')
  }

  return(variant_calls)
}
