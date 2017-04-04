.importConservation <- function(fName,
                                outName,
                                chromosome,
                                scoringType,
                                silent = FALSE,
                                writeLog = TRUE) {
    # Validate fName if file and stuff

    # variable names subject to changes

    # Setup
    open_file(fName)
    readline(n = 1)
    currPos <- extract from read line
    # Validate if the chromosome and extractedChromosomeSource are the same
    # "1" == "1"
    # "Y" == "Y"
    # "M" == "M"
    extractedChromosomeSource <- extract from read line
    # or could specify from function signature
    outName <- "chromosome + chromosomenum.RDS"

    # Whether to store coordinate with score
    chromosomeScores <- matrix size 150 million

    currSize <- size(chromosomeScores)

    while (length(entry <- readLines(opened_file, n = 1)) > 0 {
        # decide on order later
        # will there be empty data in between indices
        # refer to email
        ## No need to grow out the matrix as we will be initializing a really big matrix
        #if currIndex > currSize {
        #    chromosomeScores
        #    update currSize
        #}
        # treat as header
        if nchar(entry) > 10 {
            extract position and update currPos
        } else {
            # would it be faster to multiply then truncate at the end or along the way?
            num <- as.numeric(entry) * 1000
            # Depending on what to store
            chromosomeScores[currPos, 1] <- currPos
            chromosomeScores[currPos, 2] <- num
            currPos <- currPos + 1
        }

    }
    # remove extra indices
    # setup metadata extractedChromosomeSource, scoring type, version, data type
    save chromosomeScores as RDS with outName
    # I thought popping off the stack will lead to chromosomeScores to be garbage collected
    rm(chromosomeScores)
}
