#' Get Flagged Alleles
#'
#' Flag alleles that are present in too many samples at high variant allele frequencies as potential errors.
#'
#' @param files List of paths to all variant files from varscan pileup2cns output
#' @param exclude_cosmic_mutations Logical indicating whether or not to exclude cosmic mutations from flagged SNPs
#' @param cosmic_mutations \code{VRanges} object with position and substitution of excluded cosmic mutations
#' @param cosmic_mut_frequency Mutations with this frequency or above in the cosmic database will be excluded
#' @param memory_saving Logical. Option to save memory if you have too many samples, but takes twice as long
# @importClassesFrom VariantAnnotation VRanges SimpleVRangesList
# @importMethodsFrom S4Vectors mcols queryHits
# @importMethodsFrom GenomicRanges findOverlaps
# @importMethodsFrom BiocGenerics match
#' @export
#' @examples
#' \dontrun{
#' files <- list.files(path = "./", pattern = "sample")
#' heme_COSMIC <- load_cosmic_mutations(cosmic_mutations_path = "./heme_COSMIC.csv")
#'
#' flagged_alleles <- get_flagged_alleles(files, exclude_cosmic_mutations = TRUE,
#' cosmic_mutations = heme_COSMIC, cosmic_mut_frequency = 3)
#' }
#' @return This function returns a \code{VRanges} object with the following information:
#' \itemize{
#'	\item seqnames
#' 	\item ranges
#'	\item ref
#'	\item alt
#'	\item refdepth
#'	\item altdepth
#'	\item sampleNames
#'	\item metadata (optional)
#'	}

### TO DO LIST:
#	- Check input to ensure it is suitable
#	- Apply basic quality filter to input prior to flagging variants
# 	- Parameter to tune quantile starting point
#	- Parameter to tune fisher test cut-off

# may have to input functions as arguments
# this is kind of inefficient as we need to loop through loading the samples twice,

get_flagged_alleles <-
function(files, exclude_cosmic_mutations = FALSE, cosmic_mutations, cosmic_mut_frequency = 3, memory_saving = FALSE){

  alleles <- VariantAnnotation::VRanges()

  if(memory_saving == FALSE){

    # iterate through files to build VRanges with all alleles and VAFs
    for(samp_path in files){
      # samp name
      samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59)

      # add any extra alleles from this sample to the alleles VRanges object
      samp_alleles <- VariantAnnotation::VRanges(seqnames = GenomicRanges::seqnames(samp),
                                                 ranges = GenomicRanges::ranges(samp),
                                                 ref = VariantAnnotation::ref(samp),
                                                 alt = VariantAnnotation::alt(samp))
      alleles <- unique(append(alleles, samp_alleles))

      # add VAFs for this sample to alleles metadata
      S4Vectors::mcols(alleles)[BiocGenerics::match(samp, alleles), samp_name] = samp$VAF
      rm(samp, samp_alleles)
    }

    ### TAG THE FREQUENT SNPS
    flagged_alleles <- flag_alleles(alleles)

  } else if(memory_saving == TRUE){

    VAFs <- c()

    # iterate through files to build VRanges with all alleles
    for(samp_path in files){
      # samp name
      samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59)

      # add any extra alleles from this sample to the alleles VRanges object
      samp_alleles <- VariantAnnotation::VRanges(seqnames = GenomicRanges::seqnames(samp),
                                                 ranges = GenomicRanges::ranges(samp),
                                                 ref = VariantAnnotation::ref(samp),
                                                 alt = VariantAnnotation::alt(samp))
      alleles <- unique(append(alleles, samp_alleles))

      # store vafs
      VAFs <- append(VAFs, samp$VAF)
      rm(samp, samp_name, samp_alleles)
    }

    ### now start flagging the alleles
    VAFlen=length(VAFs)   # length of all VAFs

    Q=stats::quantile(VAFs, seq(0.99,1,0.001))
    Q=Q[-which(Q==1)]   # remove 100th percentile
    rm(VAFs)

    # matrix to store occurrences for each quantile
    vafquantile <- matrix(data = 0, nrow = length(alleles), ncol = length(Q),
                          dimnames = list(NULL, paste0("Q", gsub(x = names(Q), pattern = "%", replacement = ""))))

    for(samp_path in files){
      # samp name
      samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff_ref = 59, MAPQ_cutoff_alt = 59)

      # iterate through quantiles and tally up
      for(j in 1:length(Q)){
        ### load sample
        samp_quantile <- samp[which(samp$VAF > Q[j])]
        if(length(samp_quantile) >= 1){
          vafquantile[BiocGenerics::match(samp_quantile, alleles),j] <- vafquantile[BiocGenerics::match(samp_quantile, alleles),j] + 1
        }
      }
      rm(samp, samp_name, samp_quantile)
    }

    flag_index <- c()

    # now iterate through and determine n for fisher
    for(j in 1:length(Q)){
      n = 0; x = 1
      while (x > 0.05) {
        n = n + 1
        x = stats::fisher.test(matrix(c(n, length(files), sum(vafquantile[,j]), VAFlen), ncol = 2), conf.int = TRUE, conf.level = 0.95)[[1]]
      }
      flag_index = append(flag_index, which(vafquantile[,j] > n))
    }
    flagged_alleles <- alleles[unique(flag_index),]
  }

  # exclude the cosmic mutations
  if(exclude_cosmic_mutations == TRUE){
    # filter for frequency above 3
    cosmic_mutations <- cosmic_mutations[cosmic_mutations$hemCOSMIC_DC >= cosmic_mut_frequency]
    # removed flagged alleles overlapping with cosmic mutations
    flagged_alleles <- setdiff_VRanges(flagged_alleles, cosmic_mutations)
  }

  return(flagged_alleles)
}
