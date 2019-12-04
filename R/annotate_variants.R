#' Annotate Variants
#'
#' Wrapper function for variant annotation from CellBase through cellbaseR
#' Returns a dataframe with mutations annotated by CellBase - requires internet connection.
#'
#' @param variant_calls \code{VRanges} object from call_all_variants output
#' @param genome Genome of input object, accepts either hg19 or hg38
#' @param batchsize Size of variant batches to send to CellBase - the API will time out if processing >100 variants at a time. If it does time out, reduce the batchsize.
#' @export
#' @examples
#' \dontrun{
#' annotated_variants <- annotate_variants(variant_calls, genome="hg19")
#' }
#' @return This function returns a \code{DataFrame} object with annotations for each variant.

annotate_variants <-
function(variant_calls, genome = c("hg19", "hg38"), batchsize=80){

  # input checks
  if(length(genome) > 1){ stop("Need to specify genome, either hg19 or hg38.") }
  if(class(variant_calls) != "VRanges"){ stop('Input "variant_calls" needs to be VRanges object') }
  if(class(batchsize) != "numeric"){ stop('Input "batchsize" needs to be numeric') }
  genome = ifelse(tolower(genome) %in% c("hg19", "grch37"), "GRCh37",
                  ifelse(tolower(genome) %in% c("hg38", "grch38"), "GRCh38", genome))

  # set up variant list
  chr <- as.character(seqnames(variant_calls)) %>% stringr::str_replace("chr", "")
  pos <- start(ranges(variant_calls))
  sub <- paste(ref(variant_calls), alt(variant_calls),sep = ":")
  var_ids <- paste(chr, pos, sub, sep=":")

  # split into chunks (cellbaseR won't process >100 vars at a time)
  var_ids <- split(var_ids, ceiling(seq_along(var_ids)/batchsize))

  # set up cellbaseR with first batch
  cb <- cellbaseR::CellBaseR()
  cbparam <- cellbaseR::CellBaseParam(assembly = genome)
  annotated <- cellbaseR::getVariant(object=cb, ids=var_ids[[1]], resource="annotation", param=cbparam)

  # iterate through subsequent batches of variants and append annotated output
  if(length(var_ids) > 1){
    for(vars in var_ids[2:length(var_ids)]){
      res = cellbaseR::getVariant(object=cb, ids=vars, resource="annotation", param=cbparam)
      annotated = vctrs::vec_rbind(annotated, res)
    }
  }

  return(annotated)
}
