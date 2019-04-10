#' Generate All Models
#'
#' Generate all error models for each trinucleotide variant context. Fits error distribution with Exp model for <10,000 sequencing depth and Weibull model for >=10,000 depth.
#'
#' @param sample \code{VRanges} object of the varscan pileup2cns output annotated with variant context
#' @importFrom foreach "%dopar%"
#' @export
#' @examples
#' \dontrun{
#'
#' # Get flagged alleles and cosmic mutations
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "./heme_COSMIC.csv")
#' flagged_alleles <- get_flagged_alleles(all_samplele_names, all_samplele_paths, exclude_cosmic_mutations = TRUE,
#'     cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
#'
#' # Load and annotate samplele
#' sample <- load_as_VRanges(sample_name = "pt123",
#'     samplele_path = "./patient_123_pileup2cns", genome = "hg19", metadata = TRUE)
#' sample <- sequence_context(sample)
#' library(MafDb.gnomADex.r2.1.hs37d5)
#' annotated_sample <- annotate_MAF(varscan_output = variants,
#'     MAF_database = MafDb.gnomADex.r2.1.hs37d5, genome = "hg19")
#'
#' # Filter model input
#' sample_model_input <- filter_model_input(model_input = annotated_sample, flagged_alleles = flagged_alleles,
#'     filter_cosmic_mutations = TRUE, cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 10)
#'
#' # Generate the error models for this samplele
#' sample_models <- generate_all_models(sample = sample_model_input, plots = FALSE)
#'
#' }
#' @return This function returns a \code{dataframe} with the following information:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used (Exp or Weibull)
#'	\item Model parameter 1
#'	\item Model parameter 2
#'	}

generate_all_models <-
function(sample) {

  # Check that input is VRanges
  if(class(sample) != "VRanges"){ stop('Input "sample" needs to be VRanges object') }

  # create empty dataframe to store each flanking sequence group, the model used, and the estimated model parameters
  groups <- data.frame(FlankingSeqGroup=unique(sample$FlankingSeqGroup),model=NA,estimate1=NA,estimate2=NA)

  # Generate models for each flanking sequence group
  m <- foreach::foreach(i=1:dim(groups)[1], .combine='c', .packages = c("VariantAnnotation","fitdistrplus")) %dopar% {
    generate_model(i,sample,groups)
  }

  groups[,2:4] <- data.table::as.data.table(matrix(unlist(strsplit(m,split = " ")), ncol = 4, byrow = TRUE))[,2:4]
  groups

  for(j in 3:4){
    groups[,j] <- suppressWarnings(as.numeric(groups[,j]))
  }
  all_models <- groups

  return(all_models)
}
