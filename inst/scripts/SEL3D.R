#Pseudocode may look like actual code...please ignore...
getPdbIds <- function(targetSeq, startPos = 1 , endPos = nchar(targetSeq), convert = TRUE, silent = TRUE){
  
  #get packages

  library(bio3d)
  library(RCurl)
  library(XML)
  library(seqinr)
  
  #convert sequence into AA 
  if(convert){
    #translate needs a vector of single chars, not a string
    AAseq <- paste(translate(s2c(targetSeq)), collapse='')
    endPos <- nchar(AAseq)
  }else{ AAseq <- targetSeq }
 
  blastOutput <- blast.pdb(AAseq)

  
  #First we drop any match with ID < 40 to save runtime since we 
  #don't want to include any of these results anyway (base case)
  blastOutput <- blastOutput$hit.tbl
  numEntriesBeforeIdentityDrop <- nrow(blastOutput)
  blastOutput <- blastOutput[blastOutput$identity > 40, ]
  numEntriesAfterIdentityDrop <- nrow(blastOutput)
  
  message(sprintf("dropped %d entries because of identity <= 40", numEntriesBeforeIdentityDrop - numEntriesAfterIdentityDrop))
  
  if(nrow(blastOutput) == 0){
    return(NULL)
  }
  
  blastOutput$tax <- NA 
  blastOutput$info <- NA
  
  finalOutput <- data.frame()
  finalOutputRange <- c()
  #Pointer to keep track of position to insert alignment ranges into final output
  finalOutputRangePointer <- 1
  
  for(entry in 1:nrow(blastOutput)){
    #We're looking up taxonomy and the description of the pdbID
    #from a different server which uses slightly different formatting
    #Reformat PDB ID as output by BLAST to separate ID from chain with a "."
    RCSBid <- gsub("_", ".", blastOutput$pdb.id[entry], ignore.case = TRUE)
    entryXMLTable <- xmlParse(getURL(paste("http://www.rcsb.org/pdb/rest/describeMol?structureId=", RCSBid, sep = "")))
    
    #get taxonomy 
    
    #We don't want to drop all the non-homo sapien matches because of
    #the edge case where a match could have 100% amino acid ID for a 
    #non homo sapien which we want to include in the table.
    #So we need to get all tax information as well as the protein
    #information, without dropping anything solely based on tax
    
    tax <- getNodeSet(entryXMLTable, '//Taxonomy/@name')[[1]][[1]]
    info <- getNodeSet(entryXMLTable, '//polymerDescription/@description')[[1]][[1]]
    
    if (length(tax) > 0) {
      blastOutput$tax[entry] <- tax
    } else {
      blastOutput$tax[entry] <- NA
    }
    
    if(length(info) > 0) {
      blastOutput$info[entry] <- info
    } else {
      blastOutput$info[entry] <- NA
    }

    
    if(!intersectingPositions(finalOutputRange, blastOutput$q.start[entry], blastOutput$q.end[entry])){
    # Finally, we then find the best match that is either Homo sapiens or 100% ID.
      if((!is.na(blastOutput$tax[entry]) && blastOutput$tax[entry] == "Homo sapiens") || blastOutput$identity[entry] == 100){
        finalOutput <- rbind(finalOutput, blastOutput[entry,])
        finalOutputRange[[finalOutputRangePointer]] <- c(blastOutput$q.start[entry]:blastOutput$q.end[entry])
        finalOutputRangePointer <- finalOutputRangePointer + 1
      }
    }
  }
  
  #no matches 
  if(nrow(finalOutput) == 0){ 
    message(sprintf("all %d matches were dropped, returning NULL", nrow(blastOutput)))  
    return(NULL) 
  }
  
  finalOutput$q.start <- finalOutput$q.start + startPos - 1
  finalOutput$q.end <- finalOutput$q.end + startPos - 1
  
  for(x in 1:length(finalOutputRange)){
    finalOutputRange[[x]] <- finalOutputRange[[x]] + startPos - 1
  }
  
  minSide <- min(finalOutputRange[[1]])
  maxSide <- max(finalOutputRange[[length(finalOutputRange)]])
  
  if(endPos - maxSide > 50){
    message(sprintf("running substr from pos %d to %d", maxSide, endPos))
    maxFrame <- getPdbIds(substr(AAseq, maxSide - startPos, endPos - startPos), 
                             startPos = maxSide, endPos = endPos, convert = FALSE, silent = silent)
    finalOutput <- rbind(finalOutput, maxFrame)
  }
  
  if(minSide > 50){
    message(sprintf("running substr from pos %d to %d", startPos, minSide))
    minFrame <- getPdbIds(substr(AAseq, 1, minSide - startPos), 
                             startPos = startPos, endPos = minSide, convert = FALSE, silent = silent)
    finalOutput <- rbind(finalOutput, minFrame)
  }
  
  for(y in 1:(length(finalOutputRange) - 1)){
    leftSide <- max(finalOutputRange[[y]])
    rightSide <- min(finalOutputRange[[y+1]])
    
    if(rightSide - leftSide > 50){
      message(sprintf("running substr from pos %d to %d", leftSide, rightSide))
      interFrame <- getPdbIds(substr(AAseq, leftSide - startPos, rightSide - startPos), 
                               startPos = leftSide, endPos = rightSide, convert = FALSE, silent = silent)
      finalOutput <- rbind(finalOutput, interFrame)
    }
  }
  
  return(finalOutput)
}
#We need to store ranges of the best alignments 
#(as start (st) and end (ed) positions) so we can
#add that match to the final output without adding
#segments which overlap in the final format.
intersectingPositions <- function(r, st, ed){
  for(x in 1:length(r)){
    if(length(intersect(r[[x]], c(st:ed))) > 0){
      return(TRUE)
    }
  }
  return(FALSE)
}
# [END]
