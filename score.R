# score.R

#' Produce a gH object: a hash of genes and their associated scores
#' (likelihood that observed mutations exceed background mutations at FDR of 1%)
#'

score <- function(gX, gH, silent = FALSE, noLog = FALSE) {
    
    check parameters with .checkArgs()
    
    open gX file
    
    # for each gene, check if it already in the gH hash and if yes, do not calculate its score (in case have a repeat of a gene in gX, would just use first calculated score)
    
    for gene in gX {
        if gene already in gH {
            break
        } else {
    
            # define background mutation rate for gene
            back_mut = (missense + nonsense)/silent   # would be more accurate to use the background mutation data from Samocha et al., 2013 (http://www.nature.com/ng/journal/v46/n9/full/ng.3050.html)
            
            # calculate observed mutation rate for gene
            observed_mut = missense + nonsense
            
            
            # assume Poisson distribution (appropriate distribution, since mutation is rare event, so both the background and observed mutation rates are expected to be low)
            # determine cut-off for lambda with FDR = 0.01
            
            lower_threshold <- qpois(0.01, back_mut, lower.tail = TRUE)
            upper_threshold <- qpois(0.01, observed_mut, lower.tail = FALSE)
            
            step <- upper_threshold - lower_threshold
            
            cut_off <- upper_threshold # initialize cut-off to maximum value (upper threshold)
            for (i in seq(start = upper_threshold, end = lower_threshold, by = step/1000)){ # try 100000 possible cut-offs; can raise this even higher if we want more accuracy for FDR
                # check if resulting FDR is out of range (more than 0.01)
                if (ppois(back_mut, cut_off, lower_tail = FALSE)/ppois(observed_mut, cut_off, lower_tail = FALSE) > 0.01){
                    break
                } else {
                    cut_off <- i # if the FDR is still below 0.01, set the cut-off lower to increase sensitivity (without introducing too many more false positives)
                }
            }
            
            # determine likelihood that observed mutations exceed background mutations at an FDR of 0.01 (using the cut-off calculated above)
            
            score <- (dpois(observed_mut, cut_off)/dpois(back_mut, cut_off))
            
            # store gene and associated score in gH file
            
            gH[gene] <- score
        }
    
    }
    

}