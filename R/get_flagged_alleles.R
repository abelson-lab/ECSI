#' Get Flagged Alleles
#'
#' Flag alleles that are present in too many samples at high variant allele frequencies as potential errors.
#'
#' @param sample_names Character vector with the names of the samples
#' @param sample_paths Character vector with the paths of the samples
#' @param recurrent_mutations \code{VRanges} object with chr, pos ref, alt of frequently mutated alleles to remove from model input. \code{GRanges} objects are also accepted, in which case filtering will occur by position.
#' @param memory_saving Logical. Option to save memory if you have a lot of samples (e.g. >500 with a 16Gb RAM machine), but takes twice as long
#' @param starting_percentile Lower VAF percentile to start looking for alleles to flag. Default is 99, but can use 95 if you want to flag more alleles (more conservative)
#' @param interval VAF interval to iterate through for flagging alleles. Default is 0.001
#' @param MAPQcutoff Minimum acceptable MAPQ score; positions below this cutoff will be excluded. Default is 59
#' @export
#' @examples
#' \dontrun{
#' # get list of file names
#' file_names <- list.files(path = "./data/", pattern = "sample")
#' hemeCOSMIC_3 <- load_recurrent_mutations("example_data/COSMIC_heme_freq3.txt", genome = "hg19")
#'
#' # sample names are first 10 characters of file name
#' all_sample_names <- substr(file_names, 1, 10)
#'
#' # file paths are dir/file_name
#' all_sample_paths <- paste0("./data/", file_names)
#'
#' # get flagged alleles
#' flagged_alleles <- get_flagged_alleles(all_sample_names, all_sample_paths,
#'     recurrent_mutations = hemeCOSMIC_3, memory_saving = FALSE)
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

get_flagged_alleles <-
function(sample_names, sample_paths, recurrent_mutations = NA, memory_saving = FALSE, starting_percentile = 99, interval = 0.001, MAPQcutoff = 59){

  alleles <- VariantAnnotation::VRanges()

  if(memory_saving == FALSE){

    # iterate through files to build VRanges with all alleles and VAFs
    for(i in 1:length(sample_paths)){

      # samp path
      samp_path <- sample_paths[i]
      samp_name <- sample_names[i]

      # samp name
      #samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff = MAPQcutoff)

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
    flagged_alleles <- flag_alleles(alleles, starting_percentile = starting_percentile, interval = interval)

    ## Very different approach if memory_saving == TRUE
  } else if(memory_saving == TRUE){

    VAFs <- c()

    # iterate through files to build VRanges with all alleles
    for(i in 1:length(sample_paths)){

      # samp path
      samp_path <- sample_paths[i]
      samp_name <- sample_names[i]

      # samp name
      #samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff = MAPQcutoff)

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

    Q=stats::quantile(VAFs, seq(starting_percentile/100,1,interval))
    Q=Q[-which(Q==1)]   # remove 100th percentile
    rm(VAFs)

    # matrix to store occurrences for each quantile
    vafquantile <- matrix(data = 0, nrow = length(alleles), ncol = length(Q),
                          dimnames = list(NULL, paste0("Q", gsub(x = names(Q), pattern = "%", replacement = ""))))

    for(i in 1:length(sample_paths)){

      # samp path
      samp_path <- sample_paths[i]
      samp_name <- sample_names[i]

      # samp name
      #samp_name <- substr(samp_path,4,12)
      print(samp_name)

      # get sample as VRanges and annotate with sequence context and MAF
      samp <- load_as_VRanges(samp_name, samp_path, metadata = TRUE) %>%
        filter_MAPQ(., MAPQ_cutoff = 59)

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

  # If user provides recurrent mutations
  if(class(recurrent_mutations) %in% c("VRanges", "GRanges")){
    # remove variants that overlap with cosmic mutations
    flagged_alleles <- subtract_VRanges(flagged_alleles, recurrent_mutations)
  }

  return(flagged_alleles)
}
