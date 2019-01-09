#' Generate All Models
#'
#' Generate all error models for each trinucleotide variant context. Fits error distribution with Exp model for <10,000 sequencing depth and Weibull model for >=10,000 depth.
#'
#' @param variants \code{VRanges} object of the varscan pileup2cns output annotated with variant context
#' @param plots Logical. Whether or not to output plots of each error distribution with the model fits
#' @import doParallel
#' @import fitdistrplus
#' @import VariantAnnotation
#' @import data.table
#' @export
#' @examples
#' model_input <- filter_model_input(model_input = samp, flagged_alleles = flagged_alleles, filter_cosmic_mutations = TRUE, cosmic_mutations = hemeCOSMIC, cosmic_mut_frequency = 10)
#' samp_models <- generate_all_models(samp = model_input, plots = FALSE)
#' @return This function returns a \code{dataframe} with the following information:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used (Exp or Weibull)
#'	\item Model parameter 1
#'	\item Model parameter 2
#'	}


### TO DO LIST:
#	- Check input to ensure it is suitable
# - Is data.table even needed

# may have to input functions as arguments

generate_all_models <-
function(samp, plots = FALSE) {
  # create empty dataframe called groups for each flanking sequence group
    # include column for model (exp vs weibull)
    # column for parameter estimate 1 and 2 (values in 1 if exp, and 1 and 2 if weibull)
  groups <- data.frame(FlankingSeqGroup=unique(samp$FlankingSeqGroup),model=NA,estimate1=NA,estimate2=NA)

  # For each i (row) in groups (total 192 rows)
    # I don't know what .combine and .packages do. Look at documentation for foreach
  m <- foreach(i=1:dim(groups)[1], .combine='c', .packages = c("VariantAnnotation","fitdistrplus")) %dopar% {
    # call generate all models function with:
      # i (number of rows)
      # samp (providing data on alt alleles)
      # first column of groups (each flankingseqgroup)
    generate_model(i,samp,groups[,1], plots)
  }

  groups[,2:4] <- as.data.table(matrix(unlist(strsplit(m,split = " ")), ncol = 4, byrow = TRUE))[,2:4]
  groups

  for(j in 3:4){
    groups[,j] <- suppressWarnings(as.numeric(groups[,j]))
  }
  all_models <- groups

  return(all_models)
}
