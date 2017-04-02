accessScores <- function(fName, start, end) {
    size <- end - start + 1
    out <- numeric(length = size)
    curr = 1

    conData<- file(fName, open = "r")

    tmp <- readLines(conData, n = 1)
    header <- unlist(strsplit(tmp, " "))

    # Locate start
    # Find better later
    for (i in 1:length(header)) {
        if (startsWith(header[i], "start")) {
            startCol <- i
            break
        }
    }

    currCoordinate <- as.integer(unlist(strsplit(header[startCol], "="))[2])
    while (length(entry <- readLines(conData, n = 1)) > 0) {
        if (currCoordinate > end) {
            break
        }
        # Decide order in import implementation later
        entrySplit <- unlist(strsplit(entry, " "))
        if (entrySplit[1] == "fixedStep") {
            currCoordinate <- as.integer(unlist(strsplit(header[startCol], "="))[2])
            next
        } else {
            # turn into number
            # maybe try catch?
            num <- as.numeric(entrySplit)
            # In actual implementation, won't need this if block as we'll store everything
            if (currCoordinate >= start) {
                out[curr] <- num
                curr <- curr + 1
            }
        }
        currCoordinate <- currCoordinate + 1
    }
    close(conData)
    return(out)
}

# Can verify by hand
a <- accessScores("chrM.phastCons100way.wigFix", 5, 6)
length(a)
b <- accessScores("chrM.phastCons100way.wigFix", 16559, 16569)
length(b)
# Maybe not
c <-accessScores("chrM.phastCons100way.wigFix", 5000, 6000)
length(c)
