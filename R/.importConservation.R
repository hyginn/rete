.importConservation <- function(fName) {
    # Validate fName if file and stuff

    # variable names subject to changes

    # Setup
    open_file(fName)
    readline(n = 1)
    currPos <- extract from read line
    chromosomeNum <- extract from read line
    # or could specify from function signature
    outName <- "chromosome + chromosomenum.RDS"

    # Whether to store coordinate with score
    chromosomeScores <- vector or matrix size 20000

    currSize <- size(chromosomeScores)

    currIndex <- 1

    while (length(entry <- readLines(opened_file, n = 1)) > 0 {
        # decide on order later
        # will there be empty data in between indices
        # refer to email
        if currIndex > currSize {
            chromosomeScores
            update currSize
        }
        # treat as header
        if nchar(entry) > 10 {
            extract position and update currPos
        } else {
            # would it be faster to multiply then truncate at the end or along the way?
            num <- as.numeric(entry) * 1000
            # Depending on what to store
            chromosomeScores[currIndex, 1] <- currPos
            chromosomeScores[currIndex, 1] <- num
            currIndex <- currIndex + 1
        }

    }
    save chromosomeScores as RDS with outName
    # I thought popping off the stack will lead to chromosomeScores to be garbage collected
    rm(chromosomeScores)
}
