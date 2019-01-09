call_all_variants <-
function(samp, samp_models, call_variants = call_variants) {
  
  # Create model and model_Pvalue columns
  samp$model <- NA
  samp$model_Pvalue <- NA
  
  # Get 192 unique trinucleotide combos
  flanking_seqs <- unique(samp$FlankingSeqGroup)
  
  # for each signature in trinucleotide combo
    # compare alternate allele count with 
  m <- foreach(i=1:length(flanking_seqs), .combine='c') %dopar% {
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
  
  return(samp)
}
