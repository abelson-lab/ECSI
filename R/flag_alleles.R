#' Flag Alleles
#'
#' Flag alleles that are present in too many samples at high variant allele frequencies as potential errors.
#'
#' @param variants \code{VRanges} object from get_flagged_alleles
#' @param metadata Logical. Determine whether or not to keep the metadata (annotations) when returning flagged alleles
#' @param starting_percentile Lower VAF percentile to start looking for alleles to flag. Default is 99, but can use 95 if you want to flag more alleles (more conservative)
#' @param interval VAF interval to iterate through for flagging alleles. Default is 0.001
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

flag_alleles <-
function(variants, metadata = FALSE, starting_percentile = 99, interval = 0.001){

  # only the VAFs without chrom / pos / ref / vaf / flanking
  vars=as.matrix(S4Vectors::mcols(variants))

  # remove NAs from VAF
  vars_notNA=vars[!is.na(vars)]
  varlen=length(vars_notNA)

  # top 99th quantile (unless set by user)
  Q=stats::quantile(as.numeric(vars_notNA), seq(starting_percentile/100,1,interval))
  Q=Q[-which(Q==1)]   # remove 100th percentile

  line=c()

  # iterate through the quantiles
  for(i in 1:length(Q)) {
    # show cutoff
    print(Q[i])
    # get quantiles below this cutoff
    q = as.data.frame(Q[which(Q <= Q[i])])
    # save cutoff percentile
    z = rownames(q)[dim(q)[1]]
    # save cutoff numeric value
    q = as.numeric(q[which(rownames(q) == z), ])

    # index of VAFs above that quantile
    index <- which(vars > q, arr.ind = T)

    # summarise index
    tab = as.data.frame(table(index[,1]))
    tab$Var1=as.numeric(as.character(tab$Var1))   # make sure Var1 is numeric
    n = 0
    x = 1

    while (x > 0.05) {
      n = n + 1
      print(n)
      x = stats::fisher.test(matrix(c(n, dim(vars)[2], sum(tab$Freq), varlen), ncol = 2), conf.int = TRUE, conf.level = 0.95)[[1]]
    }

    print(n)
    index = which(tab$Freq > n)
    print(length(index))
    line = append(line, as.numeric(as.character(tab$Var1[index])))
  }

  # report flagged alleles
  flagged_alleles <- variants[unique(line),]

  ## Remove the metadata to save memory
  if (metadata == FALSE) {
    S4Vectors::mcols(flagged_alleles) <- NULL
  }

  return(flagged_alleles)
}
