#' Generate Model
#'
#' Generate error model given a specific trinucleotide context. Fits error distribution with Exp model for <10,000 sequencing depth and Weibull model for >=10,000 depth.
#'
#' @param i Index of the FlankingSeqGroup for which the model is being generated
#' @param data \code{Dataframe} from generate_all_models with FlankingSeqGroup and three empty columns for parameters of model fit
#' @param groups \code{Dataframe} with 192 rows and 4 columns. First column contains the 192 trinucleotide contexts/FlankingSeqGroup
#' @param plot Logical. Whether or not to output a plot of the error distribution with the model fit
# @importFrom fitdistrplus fitdist
# @importFrom graphics hist title points
# @importClassesFrom VariantAnnotation VRanges
# @importMethodsFrom VariantAnnotation altDepth refDepth
#' @return This function returns a \code{dataframe} with the following information:
#' \itemize{
#'	\item FlankingSeqGroup
#' 	\item Model Used (Exp or Weibull)
#'	\item Model parameter 1
#'	\item Model parameter 2
#'	}

### TO DO LIST:
#	- Check input to ensure it is suitable
#	- Reduce the number of comments
#	- Ammend the generated plot

generate_model <-
function(i, data, groups, plot=FALSE){

  # set up output for function.
    # i, fit, param estimate 1, param estimate 2
  l <- paste(c(i,"None",NA,NA))

  # identify the FlankingSeqGroup in the data corresponding to index i in groups
  index <- which(data$FlankingSeqGroup == groups$FlankingSeqGroup[i])

  # tabulate alternate allele counts for this group into contingency table
  tab <- as.data.frame(table(VariantAnnotation::altDepth(data)[index]))
  # change alternate allele counts to numeric
    ### TO DO: Is this necessary???
  tab$Var1 <- as.numeric(as.character(tab$Var1))

  # get alternate allele count
  altBases <- VariantAnnotation::altDepth(data)[index]

  # Display Alt Allele count histogram and log(VAF) ~ log(Coverage)
  if(plot==TRUE){
    # histogram of alt allele count
    graphics::hist(VariantAnnotation::altDepth(data)[index],breaks = 100)
    # contingency table of alt allele count
    table(altBases)

    # get data corresponding to this signature
    s <- data[index,]
    # establish plot x and y limits
    xmin <- min(log(VariantAnnotation::refDepth(s) + VariantAnnotation::altDepth(s)))
    xmax <- max(log(VariantAnnotation::refDepth(s) + VariantAnnotation::altDepth(s)))
    ymin <- min(log(VariantAnnotation::altDepth(s) / (VariantAnnotation::refDepth(s) + VariantAnnotation::altDepth(s))))
    ymax <- max(log(VariantAnnotation::altDepth(s) / (VariantAnnotation::refDepth(s) + VariantAnnotation::altDepth(s))))

    # Plot the relationship between log(VAF) ~ log(Depth) for that signature
    plot(x = log(VariantAnnotation::refDepth(s)+VariantAnnotation::altDepth(s)),
         y = log(VariantAnnotation::altDepth(s)/(VariantAnnotation::refDepth(s)+VariantAnnotation::altDepth(s))),
         xlab = "log( Read Depth )", ylab = "log( Variant Allele Frequency )", xlim = c(xmin,xmax), ylim = c(ymin,ymax))
    graphics::title(paste("log(VAF) ~ log(Depth) for Signature:", groups$FlankingSeqGroup[i]))

    # colour in point with max Alternate Allele Frequency
      ### TO DO: Is this necessary???
    max_alt_count <- which(VariantAnnotation::altDepth(s) == max(VariantAnnotation::altDepth(s)))
    graphics::points(x = log(VariantAnnotation::refDepth(s)+VariantAnnotation::altDepth(s))[max_alt_count],
           y = log(VariantAnnotation::altDepth(s)/(VariantAnnotation::refDepth(s)+VariantAnnotation::altDepth(s)))[max_alt_count], pch=19, col="red")
    }

  # If there are more than two different alternate allele counts for this seq group
  if(length(unique(altBases))>2){

    # divide each alt allele count by immediately smaller alt count
    index1 <- tab$Var1[-1]/tab$Var1[-length(tab$Var1)]     # this corresponds to the ratios

    # which of these ratios is greater than 1.41
    index2 <- which(index1>1.41)     # index of the ratios > 1.41

    # get the numerator from each ratio > 1.41
    index3 <- tab$Var1[index2+1]     # numerators of the ratios > 1.41

    # what is the ratio > 1.41 with the smallest numerator from above 14?
    index4 <- min(which(index3>=14))     # index in index2 of ratio >1.41 with smallest numerator > 14

    # what is the index of that numerator in the original contingency table
    index5 <- index2[index4]+1     # index of smallest number above 14 AND with ratio > 1.41 in contingency table

    # if a number > 14 AND with a ratio > 1.41 exists
    if(!is.na(index5)){
      # remove that number and every number larger than it from altBases
      altBases <- altBases[-match(tab$Var1[c(index5:dim(tab)[1])],altBases)]
      # remove that number and every number larger than it from tab (contingency table)
      tab <- tab[-c(index5:dim(tab)[1]),]
    }

    # Temp table with alt allele count from min value to max value.
    # need to convert tab from factor -> character -> numeric
    # Frequency column is empty.
    tabTmp <- data.frame(Var1=c(min(as.numeric(as.character(tab$Var1))):max(as.numeric(as.character(tab$Var1)))),Freq=0)
    # Fill in temp table with frequencies from tab
    tabTmp$Freq <- tab$Freq[match(tabTmp$Var1,tab$Var1)]   # this is a Vlookup function for Var1
    # Replace NA with 0
    tabTmp$Freq[which(is.na(tabTmp$Freq))] <- 0

    # Replace tab with tabTmp
    tab <- tabTmp
    # save alternate allele counts as x1?
    x1 <- as.numeric(tab$Var1)

    #weibull
    fitW <- try(fitdistrplus::fitdist(altBases, distr = "weibull",method = "qme",probs=c(0.5,0.99)),silent=TRUE)
    #exponential
    fitEx <- try(fitdistrplus::fitdist(altBases, distr = "exp",method = "qme",probs=c(0.99)),silent=TRUE)

    # if mean reads less than 10,000
    if(stats::median(VariantAnnotation::refDepth(data))<=10000){
      # make sure that exponential fit worked
      if(!inherits(fitEx, "try-error") & !is.na(fitEx$estimate[[1]])){
        # Ex as empty matrix with 2 cols
        Ex <- matrix(nrow = 0, ncol = 2)
        # iterate through from 0.05 to Exp Fit rate estimate + 0.5, taking steps of 0.05
        for(r in seq(0.05,fitEx$estimate[[1]]+0.5,0.05)){
          # generate exponential distribution of alt count for each decay rate
          Density <- stats::dexp(x1, rate = r)

          # vector combining observed Alt Count with expected Alt Count from Density plot
          x <- cbind(tab$Freq,(tab$Freq[1]/Density[1])*Density)
          # ratio of Observed Alt Count over Expected Alt Count
          a <- x[,1]/x[,2]
          # find ratios less than 1: expected decay rate slower than observed
          index <- which(a<1)
          # if these ratios exist, get reciprocal of the ratios that are less than 1.
          if(length(index>0)){a[which(a<1)] <- 1/a[which(a<1)]}
          # if any of these ratios are now infinity (i.e. 0 observed alt count):
          index <- which(a=="Inf")
          # remove the ones that are infinity
          if(length(index>0)){a <- a[-index]}

          # get the sum of the ratios between expected and observed.
          # the magnitude of the ratio represents discrepancy between expected and observed
          # the best fit will be the smallest ratio (closest to 1)
          sABS <- sum(a)

          # append this to the Exponential info for this SequenceContext.
          # can use to identify optimal rate for exponential decay that fits the data
          Ex <- rbind(Ex, c(r, sABS))
        }
      }

      # Turn Ex into a data frame
      Ex <- as.data.frame(Ex)
      # Convert Ex columns (factors) into numeric
      #for(j in 1:4){
      for(j in 1:2){
        Ex[,j] <- as.numeric(as.character(Ex[,j]))
      }
      # Which discrepancy ratio is the lowest (model best fits the data)
      indexEx <- which(Ex$V2 == min(Ex$V2, na.rm = T))
      # If there is more than one ideal model:
      if(length(indexEx)>1){
        # fit the distribution and get a rate estimate
        ### Fit the quantiles by a 0.9 probability vector
        ### I DON'T KNOW WHAT THIS MEANS OR HOW TO OPTIMIZE IT ##################
        fitEx <- fitdistrplus::fitdist(altBases, distr = "exp",method = "qme",probs=c(0.9))

        ##### I DON'T UNDERSTAND THIS ONE EITHER #############################
        ### Do we just force in a rate estimate?
        ### Why not choose between the available estimates?
        ### Why not do this from the get go?
        Ex$V1[i] <- fitEx$estimate[[1]]
        indexEx <- i
      }
    }

    # Choose whether to use Exp or Weibull based on sequence depth
    if(stats::median(VariantAnnotation::refDepth(data))<=10000){
      l=paste(c(i,"exp",Ex$V1[indexEx], "NA" )) # exp for depth < 10000
    } else {
      l=paste(c(i,"weibull",fitW$estimate[[1]],fitW$estimate[[2]])) # weibull for depth > 10000
    }
  }

  # return l vector with i, model, estimate1, estimate
  return(l)
}
