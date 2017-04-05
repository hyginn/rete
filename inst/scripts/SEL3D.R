#SEL3D.R
#
#' \code{SEL3D} wrapper function for getPdbIds. Blasts input sequence and
#'              return MT object with match information.
#'
#' Details.
#'
#' @param targetSeq input sequence to be queried, must be cDNA of a length divisible by 3.
#' @param convert bool to convert DNA sequences to amino acids, True by default
#'                NOTE: The user should NOT change convert unless for testing purposes
#' @param silent bool, set to True to turn off all messages.
#' @param writeLog bool, default set to TRUE to write a log
#' @return MT object - dataframe with headers queryid, subjectids, identity, positives,
#'                     alignment length, mismatches, gapopens, q.start, q.end, s.start, s.end,
#'                     evalue, bitscore, mlog.evalue, pdb.id
#'
#'
#'
#' @examples
#' >SEL3D("AAGTCT")
#' >sequence = "ACTGTGACTGTAAATGTT"
#' >SEL3D(sequence)
#'

SEL3D <- function(targetSeq, convert = TRUE, silent = TRUE, writeLog = TRUE) {
    #verify params
    if(mode(convert) != "logical") {
        stop(paste0("Mode error: argument requires logical input"))
    }
    if(mode(silent) != "logical") {
        stop(paste0("Mode error: argument requires logical input"))
    }
    return(getPdbIds(targetSeq, convert = convert, silent = silent, writeLog = writeLog))
}

#getPdbIds
#
#' \code{getPdbIds} Blasts targetSeq to PDB and returns an MT object

#' Details.
#'
#' @param targetSeq input sequence to be queried, must be cDNA of a length divisible by 3.
#' @param startPos start position of the sequence, set to 1 by default, altered in recursive steps
#' @param endPos end position of sequence, set to 0, altered to aa length in recursive steps
#' @param convert bool to convert DNA sequences to amino acids, True by default
#' @param silent bool, set to True to turn off all messages.
#' @return MT object - dataframe with headers queryid, subjectids, identity, positives,
#'                     alignment length, mismatches, gapopens, q.start, q.end, s.start, s.end,
#'                     evalue, bitscore, mlog.evalue, pdb.id
#'
#'
#'
#' @examples
#' >getPdbIds("AAGTCT")
#' >sequence = "ACTGTGACTGTAAATGTT"
#' >getPdbIds(sequence)
#' >getPdbIds("MGEVST", convert=FALSE)

getPdbIds <- function(targetSeq, startPos = 1 , endPos = nchar(targetSeq), convert = TRUE, silent = TRUE, writeLog = TRUE){


    #check for valid targetSeq
    if(targetSeq =="") {
        stop(paste0("Error: empty targetSeq"))
    }
    if(convert && !nchar(targetSeq) %% 3 == 0) {
        stop(paste0("Error: targetSeq's length needs to be divisible by 3"))
    }
    #we want to blast sequences of amino acid length >50 so the input cDNA length should be >150
    if(nchar(targetSeq) < 150) {
        stop(paste0("Error: targetSeq's length must be >150"))
    }
    if(convert) {
        if(grepl("[^AGTC]", targetSeq)) {
            stop(paste0("Error: targetSeq contains invalid characters"))}
    }

    #convert sequence into AA
    if(convert){
        #translate needs a vector of single chars, not a string
        AAseq <- paste(seqinr::translate(seqinr::s2c(targetSeq)), collapse='')
        endPos <- nchar(AAseq)
    } else { AAseq <- targetSeq }


    #log file arguments
    call_params <- c()
    call_params[1] <- sprintf("sequence run with positions %d to %d ", startPos, endPos)
    call_params[2] <- sprintf("convert: %s ", as.character(convert))
    call_params[3] <- sprintf("silent: %s ", as.character(silent))
    parameters <- paste0(call_params, collapse = ",")

    logTitle <- "SEL3D"

    #local blast

    #we have to generate a file for blastp to work
    write(AAseq, file = "inst/extdata/tmp/sequence.txt", append = FALSE)

    #this is the run command
    local_blastp_command <- "blastp -db inst/extdata/pdbaa/pdbaa -query inst/extdata/tmp/sequence.txt -outfmt 10"

    #this is getting the output
    local_blast_output <- as.data.frame(read.csv(pipe(local_blastp_command), header = FALSE))

    #returns csv with no columns, so manual add here
    col_labels <- c("queryid",
                    "pdb.id",
                    "identity",
                    "alignmentlength",
                    "mismatches",
                    "gapopens",
                    "q.start",
                    "q.end",
                    "s.start",
                    "s.end",
                    "evalue",
                    "bitscore")
    #setting columns
    colnames(local_blast_output) <- col_labels


    #non local blast
    # if(silent) {
    #   suppressMessages(capture.output(blastOutput <- bio3d::blast.pdb(AAseq)))
    # } else {
    #   blastOutput <- bio3d::blast.pdb(AAseq)
    # }
    # blastOutput <- blastOutput$hit.tbl

    blastOutput <- local_blast_output
    #First we drop any match with ID < 40 to save runtime since we
    #don't want to include any of these results anyway (base case)

    numEntriesBeforeIdentityDrop <- nrow(blastOutput)
    blastOutput <- blastOutput[blastOutput$identity > 40, ]
    numEntriesAfterIdentityDrop <- nrow(blastOutput)

    if(!silent) {
        message(sprintf("dropped %d entries because of identity <= 40",
                        (numEntriesBeforeIdentityDrop - numEntriesAfterIdentityDrop)))
    }

    if(writeLog){
        logEvent(eventTitle = logTitle, eventCall = parameters,
                 notes = sprintf("dropped %d entries because of identity <= 40",
                                 numEntriesBeforeIdentityDrop - numEntriesAfterIdentityDrop))
    }

    if(nrow(blastOutput) == 0){
        if(writeLog){
            logEvent(eventTitle = logTitle, eventCall = parameters,
                     notes = sprintf("dropped all entries, returning NULL"))
        }
        if(!silent){
            message("dropped all entries, returning NULL")
        }
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
        entryXMLTable <- XML::xmlParse(RCurl::getURL(paste
                                                     ("http://www.rcsb.org/pdb/rest/describeMol?structureId=",
                                                         RCSBid, sep= "")))

        #get taxonomy

        #We don't want to drop all the non-homo sapien matches because of
        #the edge case where a match could have 100% amino acid ID for a
        #non homo sapien which we want to include in the table.
        #So we need to get all tax information as well as the protein
        #information, without dropping anything solely based on tax

        tax <- XML::getNodeSet(entryXMLTable, '//Taxonomy/@name')[[1]][[1]]
        info <- XML::getNodeSet(entryXMLTable, '//polymerDescription/@description')[[1]][[1]]

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
            # Finally, we then find the best match that is Homo sapiens
            if((!is.na(blastOutput$tax[entry]) && blastOutput$tax[entry] == "Homo sapiens")){
                finalOutput <- rbind(finalOutput, blastOutput[entry,])
                finalOutputRange[[finalOutputRangePointer]] <-
                    c(blastOutput$q.start[entry]:blastOutput$q.end[entry])
                finalOutputRangePointer <- finalOutputRangePointer + 1
            }
        }
    }

    #remove all alignments with any gaps
    finalOutput <- finalOutput[finalOutput$gapopens == 0,]

    #no matches were found, print how many were dropped
    if(nrow(finalOutput) == 0){
        if(!silent) {
            message(sprintf("all %d matches were dropped, returning NULL", nrow(blastOutput)))
        }
        if(writeLog){
            logEvent(eventTitle = logTitle, eventCall = parameters,
                     notes = sprintf("all %d matches were dropped, returning NULL", nrow(blastOutput)))
        }
        return(NULL)
    }

    #For the recursion, we need to store the start and end positions of our new
    #input sequence with reference to the original sequence.  However, the recursion
    #will just think it's a new sequence starting from index 1 to its length
    #so after we blast the subsequence, we need to adjust the qstart and
    #qend positions by the numbers we stored before the recursion.

    finalOutput$q.start <- finalOutput$q.start + startPos - 1
    finalOutput$q.end <- finalOutput$q.end + startPos - 1

    #Check case where first aligned sequence is in the middle of the sequence
    minSide <- 0
    maxSide <- 0

    for(x in 1:length(finalOutputRange)){
        finalOutputRange[[x]] <- finalOutputRange[[x]] + startPos - 1
        if(minSide == 0 || min(finalOutputRange[[x]] < minSide)){
            minSide <- min(finalOutputRange[[x]])
        }

        if(maxSide == 0 || max(finalOutputRange[[x]] > maxSide)){
            maxSide <- max(finalOutputRange[[x]])
        }
    }


    #Recursive step: if there are subsequences remaining after our right-most query that
    #are >50aa in length, blast these as well

    #Check for gaps between last aligned amino acid and last amino acid of
    #this recursive step's sequence
    if(endPos - maxSide > 50){
        if(!silent) {
            message(sprintf("running substr from pos %d to %d", maxSide, endPos))
        }
        maxFrame <- getPdbIds(substr(AAseq, maxSide - startPos, endPos - startPos),
                              startPos = maxSide, endPos = endPos, convert = FALSE, silent = silent)
        finalOutput <- rbind(finalOutput, maxFrame)
    }

    #Recursive step: if there are subsequences remaining before our leftmost-most query that
    #are >50aa in length, blast these as well
    if(minSide - startPos > 50){
        if(!silent) {
            message(sprintf("running substr from pos %d to %d", startPos, minSide))
        }
        minFrame <- getPdbIds(substr(AAseq, 1, minSide - startPos),
                              startPos = startPos, endPos = minSide, convert = FALSE, silent = silent)
        finalOutput <- rbind(finalOutput, minFrame)
    }

    #We've covered the recursive steps which blast squences book-ending our first and last alignments
    #now we have to check subsequences between our alignments (if there are multiple alignments) for
    #the current recursive step

    if(length(finalOutputRange) > 1){

        #If the best alignment is not the leftmost, sort by the first postion that
        #an alignment exists to look for gaps
        finalOutputRange <- finalOutputRange[order(sapply(finalOutputRange,'[[',1))]

        for(y in 1:(length(finalOutputRange) - 1)){
            leftSide <- max(finalOutputRange[[y]])
            rightSide <- min(finalOutputRange[[y+1]])

            if(rightSide - leftSide > 50){
                if(!silent) {
                    message(sprintf("running substr from pos %d to %d", leftSide, rightSide))
                }
                interFrame <- getPdbIds(substr(AAseq, leftSide - startPos, rightSide - startPos),
                                        startPos = leftSide, endPos = rightSide, convert = FALSE, silent = silent)
                finalOutput <- rbind(finalOutput, interFrame)
            }
        }
    }

    if(writeLog){
        logEvent(eventTitle = logTitle, eventCall = parameters,
                 notes = sprintf("matched %d pdb ids to target sequence", nrow(finalOutput)))
    }
    return(finalOutput)
}
#intersectingPositions
#
#' \code{intersectingPositions} given a list of ranges and a new range, if the new
#'                              range intersects with any range in the list,
#'                              return TRUE, else FALSE
#' Details.
#'
#' @param r list of ranges to check whether range st-ed intersects with
#' @param st start of new range
#' @param ed end of new range
#' @return bool
#'
#' @examples
#' > intersectingPositions(c(2,4,6), 5, 7)
#' [1] TRUE
#'
intersectingPositions <- function(r, st, ed){

    for(x in 1:length(r)){
        if(length(intersect(r[[x]], c(st:ed))) > 0){
            return(TRUE)
        }
    }
    return(FALSE)
}
# [END]
