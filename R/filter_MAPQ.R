filter_MAPQ <-
function(varscan_output, MAPQ_cutoff = 59) {
  # keep variants above MAQ cutoff
  varscan_output <- varscan_output[which(varscan_output$MapQual2 >= MAPQ_cutoff),]

  return(varscan_output)
}
