annotate_MAF <-
function(varscan_output, liftOver_chain_hg19toHg38 = "hg19ToHg38.over.chain", liftOver_chain_hg38toHg19 = "hg38ToHg19.over.chain"){
  
  # Liftover from hg19 to hg38
  hg19tohg38 <- import.chain(liftOver_chain_hg19toHg38)
  varscan_output_hg38 <- unlist(liftOver(varscan_output, hg19tohg38))
  genome(varscan_output_hg38) <- "GRCh38"

  # Annotate with MAF from Exac
  varscan_output_hg38_anno <- gscores(MafDb.gnomAD.r2.0.1.GRCh38, varscan_output_hg38)
  seqlevelsStyle(varscan_output_hg38_anno) <- "UCSC"
  
  # Liftover from hg38 back to hg19
  hg38tohg19 <- import.chain(liftOver_chain_hg38toHg19)
  varscan_output_anno <- unlist(liftOver(varscan_output_hg38_anno, hg38tohg19))
  genome(varscan_output_anno) <- "hg19"
  
  return(varscan_output_anno)
}
