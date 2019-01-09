generate_model <-
function(i, data, FlankingSeqGroup, plots=FALSE){
  
  # set up output for function. 
    # i, fit, param estimate 1, param estimate 2
  l <- paste(c(i,"None",NA,NA))
  
  # identify the FlankingSeqGroup in the data corresponding to index i in groups
  index <- which(data$FlankingSeqGroup==groups$FlankingSeqGroup[i])
  
  # tabulate alternate allele counts for this group into contingency table
  tab <- as.data.frame(table(altDepth(data)[index]))
  # change alternate allele counts to numeric
    ### TO DO: Is this necessary???
  tab$Var1 <- as.numeric(as.character(tab$Var1))
  
  # get alternate allele count
  altBases <- altDepth(data)[index]
  
  # Display Alt Allele count histogram and log(VAF) ~ log(Coverage)
  if(plots==TRUE){
    # histogram of alt allele count
    hist(altDepth(data)[index],breaks = 100)
    # contingency table of alt allele count
    table(altBases) 
    
    # get data corresponding to this signature
    s <- data[index,] 
    # establish plot x and y limits
    xmin <- min(log(refDepth(s)+altDepth(s))) 
    xmax <- max(log(refDepth(s)+altDepth(s))) 
    ymin <- min(log(altDepth(s)/(refDepth(s)+altDepth(s)))) 
    ymax <- max(log(altDepth(s)/(refDepth(s)+altDepth(s)))) 

    # Plot the relationship between log(VAF) ~ log(Depth) for that signature
    plot(x = log(refDepth(s)+altDepth(s)), y = log(altDepth(s)/(refDepth(s)+altDepth(s))),
         xlab = "log( Read Depth )", ylab = "log( Variant Allele Frequency )", xlim = c(xmin,xmax), ylim = c(ymin,ymax)) 
    title(paste("log(VAF) ~ log(Depth) for Signature:", groups$FlankingSeqGroup[i]))

    # colour in point with max Alternate Allele Frequency
      ### TO DO: Is this necessary???
    max_alt_count <- which(altDepth(s) == max(altDepth(s))) 
    points(x = log(refDepth(s)+altDepth(s))[max_alt_count], y = log(altDepth(s)/(refDepth(s)+altDepth(s)))[max_alt_count], pch=19, col="red")
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
    fitW <- try(fitdist(altBases, distr = "weibull",method = "qme",probs=c(0.5,0.99)),silent=TRUE)
    #exponential
    fitEx <- try(fitdist(altBases, distr = "exp",method = "qme",probs=c(0.99)),silent=TRUE)
    
    # if mean reads less than 10,000
    if(median(refDepth(data))<=10000){
      # make sure that exponential fit worked
      if(!inherits(fitEx, "try-error") & !is.na(fitEx$estimate[[1]])){
        # Ex as empty matrix with 4 cols 
        Ex <- matrix(nrow = 0, ncol = 4)
        # iterate through from 0.05 to Exp Fit rate estimate + 0.5, taking steps of 0.05
        for(r in seq(0.05,fitEx$estimate[[1]]+0.5,0.05)){
          # generate exponential distribution of alt count for each decay rate
          Density <- dexp(x1, rate = r)
          
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
        fitEx <- fitdist(altBases, distr = "exp",method = "qme",probs=c(0.9))
        
        ##### I DON'T UNDERSTAND THIS ONE EITHER #############################
        ### Do we just force in a rate estimate? 
        ### Why not choose between the available estimates?
        ### Why not do this from the get go?
        Ex$V1[i] <- fitEx$estimate[[1]]
        indexEx <- i
      }
    }
    
    # Choose whether to use Exp or Weibull based on sequence depth
    if(median(refDepth(data))<=10000){
      best <- 1
    } else {
      best <- 2
    }
    # if <= 10,000 use the refined Exp dist
    if(best==1){l=paste(c(i,"exp",Ex$V1[indexEx], "NA" ))}
    # weibull is super easy eh?
    if(best==2){l=paste(c(i,"weibull",fitW$estimate[[1]],fitW$estimate[[2]]))}
  }
   
  # return l vector with i, model, estimate1, estimate
  return(l)
}
