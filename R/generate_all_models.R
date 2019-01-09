generate_all_models <-
function(samp, generate_model = generate_model, plots = FALSE) {
  # create empty dataframe called groups for each flanking sequence group
    # include column for model (exp vs weibull) 
    # column for parameter estimate 1 and 2 (values in 1 if exp, and 1 and 2 if weibull)
  groups <- data.frame(FlankingSeqGroup=unique(samp$FlankingSeqGroup),model=NA,estimate1=NA,estimate2=NA)
  
  # For each i (row) in groups (total 92 rows)
    # I don't know what .combine and .packages do. Look at documentation for foreach
  m <- foreach(i=1:dim(groups)[1], .combine='c', .packages = "fitdistrplus") %dopar% {
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
