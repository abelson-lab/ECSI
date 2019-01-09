filter_model_input <-
function(model_input, flagged_alleles, MAF_cutoff = 0.001, VAF_cutoff = 0.05, 
                               filter_cosmic_mutations = FALSE, cosmic_mutations, cosmic_mut_frequency = 10) {
  # keep variants below MAF cutoff
  model_input <- model_input[which(model_input$AF < MAF_cutoff),]
  # keep variants below VAF cutoff
  model_input <- model_input[which(model_input$VAF < VAF_cutoff),]
  
  # remove variants that overlap with flagged alleles
  model_input <- model_input[-queryHits(findOverlaps(model_input, flagged_alleles, type = "any"))]
  
  # if filtering cosmic mutations 
  if(filter_cosmic_mutations == TRUE){
    # filter for frequency above 10
    cosmic_mutations <- cosmic_mutations[cosmic_mutations$hemCOSMIC_DC >= cosmic_mut_frequency]
    # remove variants that overlap with cosmic mutations
    model_input <- model_input[-queryHits(findOverlaps(model_input, cosmic_mutations, type = "any"))]
  }
  
  return(model_input)
}
