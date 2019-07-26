#' Call Variants
#'
#' Call all variants from the varscan pileup2cns output given a specifci trinucleotide context and the corresponding error model.
#'
#' @param i Index of the FlankingSeqGroup for which the model is being generated
#' @param data \code{Data.Frame} from \code{generate_all_models} with FlankingSeqGroup and three empty columns for parameters of model fit
#' @param flanking_seqs String with one of 192 trinucleotide contexts in which variants may be found
#' @param error_models \code{Data.Frame} of context specific error models generated from \code{generate_all_models}
#' @return This function returns a \code{VRanges} with the following metadata:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used
#'	\item Model pvalue
#'	}


call_variants <-
function(i, data, flanking_seqs, error_models){

  # Define empty vectors
  p <- c()
  ind1 <- c()
  ind2 <- c()
  ind3 <- c()

  # Get corresponding model for flanking seq i
  index <- which(error_models[,1] == flanking_seqs[i])

  # Get corresponding varscan calls for flanking seq i
  index1 <- which(data$FlankingSeqGroup==flanking_seqs[i])

  # Get unique alt allele counts for flanking seq i
  alt_counts <- unique(VariantAnnotation::altDepth(data)[index1])

  # If index exists
  if(flanking_seqs[i] %in% error_models[,1]){

    ### If error model exists for this flanking seq
    if(error_models$model[index] != "None"){

      # If error model for this flanking seq is exp (<10,000 mean depth)
      if (error_models$model[index] == "exp") {
        # Get estimate for exp decay rate
        estimate1 <- error_models$estimate1[index]

        # For each unique alt allele count
        for(j in alt_counts){
          # Compare alt allele count against error distribution
            # concatenate into vector with previous
          p <- c(p,1-stats::pexp(j,estimate1))
        }

      # If error model for this flanking seq is weibull (>10,000 mean depth)
      } else if (error_models$model[index] == "weibull") {
        # Get estimates for weibull parameters scale and shape
        estimate1 <- error_models$estimate1[index]
        estimate2 <- error_models$estimate2[index]

        # For each unique alt allele count
        for(j in alt_counts){
          # get pvalue from error dist for this alt allele count
          p <- c(p,1-stats::pweibull(j,estimate1,estimate2))
        }
      }

      # iterate through alt allele counts
      for(z in 1:length(alt_counts)){
        # find the rows in the varscan call corresponding to this alt allele count
        index2 <- which(VariantAnnotation::altDepth(data)[index1] == alt_counts[z])
        # index for rows in varscan call corresponding to this alt count + Flanking Seq Group
        ind1 <- append(ind1,index1[index2])
        # p_value corresponding to each entry in ind1
        ind2 <- append(ind2,rep(p[z],length(index1[index2])))
        # model used corresponding to each entry in ind1
        ind3 <- append(ind3,rep(error_models$model[index],length(index1[index2])))
      }

    ### If no error model exists for this flanking seq
      # assign NA for now. will replace with VarScan Pvalue
    } else{
      for(z in 1:length(alt_counts)){
        index2 <- which(VariantAnnotation::altDepth(data)[index1]==alt_counts[z])
        ind1 <- append(ind1,index1[index2])
        ind2 <- append(ind2,rep(NA,length(index1[index2])))
        ind3 <- append(ind3,rep("None",length(index1[index2])))
      }
    }
  # if index doesn't exist, also make pvalue NA
  } else{
    for(z in 1:length(alt_counts)){
      index2 <- which(VariantAnnotation::altDepth(data)[index1]==alt_counts[z])
      ind1 <- append(ind1,index1[index2])
      ind2 <- append(ind2,rep(NA,length(index1[index2])))
      ind3 <- append(ind3,rep("None",length(index1[index2])))
    }
  }

  return(list(ind1,ind2,ind3))
}
