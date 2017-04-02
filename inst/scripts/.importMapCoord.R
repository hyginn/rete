# importMapCoord.R
# non-exported

#==========================makeMap and Validation===============================

# Helper function - generates the closure for a single gene model
makeMap <- function(ucscMod) {

    # Check if the number of exons listed in start and end match
    if(length(ucscMod$exonStart) != length(ucscMod$exonEnd)) {
        return(NULL)
    }

    # Prepare the environment (exon matrix, start, end, strand)
    exons <- cbind(ucscMod$exonStart, ucscMod$exonEnd)
    start <- ucscMod$start
    end <- ucscMod$end
    reverse <- FALSE
    if (ucscMod$strand == "-") {
        reverse <- TRUE
    }

    # Check if the exons and start stop points are validly ordered
    if(is.unsorted(c(start, as.vector(t(exons)), end))) {
        return(NULL)
    }

    # Return a function that uses the exon matrix, strand and start, stop positions
    return(
        function(x) {
            if (x < start | x > end) {
                # out-of-range coordinate
                return(NULL)
            }

            # Determine the exon our coordinate is located in
            iEx <- which(x >= exons[ , 1] & x <= exons[ , 2])
            if(length(iEx) != 1){
                # intron
                return(0)
            }

            # Determine the number of bases in previous exons
            if (iEx == 1) {
                lPre <- 0
            } else {
                lPre <- sum((exons[ , 2] - exons [ , 1] + 1)[1:(iEx-1)])
            }

            # Determine the number of bases from the current
            lIn <- x - exons[iEx, 1] + 1

            if (reverse) {
                return(sum(exons[ , 2] - exons [ , 1] + 1) - lPre - lIn + 1)
            } else {
                return(lPre + lIn)
            }
        }
    )
}



# Generate random start, stop, intron, exon, and strand data for the gene model
# Ideally, we would use the importGeneData function to get real data
# when it gets implemented

# Convert a sorted even-sized list with a random number of random numbers
# within an interval into our "gene model"

set.seed(3141592)

lG <- 2*sample(2:10, 1)
preMod <- sample(1:10000, size=lG)
preMod <- sort(preMod)
plusmin <- c('+', '-')

# The first and last integers become the start and stop coordinates
# While the middle integers are paired up to store exon start-stop coordinates
preMod

geneMod <- list(start=preMod[1], end=preMod[lG], strand=plusmin[sample(1:2, 1)])
geneMod$exonStart <- c(preMod[2:(lG-1)])[c(TRUE, FALSE)]
geneMod$exonEnd <- c(preMod[3:(lG)])[c(TRUE, FALSE)]

sampleFx <- makeMap(geneMod)

# Regenerate the gene model from the results
expExon <- vector(mode="numeric", length=lG-2)


# Loop through all bases in the gene
previous <- 0
counter <- 1
for (i in preMod[1]:preMod[lG]) {
    current <- sampleFx(i)
    if (previous == 0 & current != 0) {
        # start of exon
        expExon[counter] <- i
        counter <- counter + 1
    }
    if (previous != 0 & current == 0) {
        # 1 after end of exon
        expExon[counter] <- i - 1
        counter <- counter + 1
    }
    previous <- current
}

expGeneMod <- c(preMod[1], expExon, preMod[lG])

# Are they the same? (They should be)
preMod
expGeneMod

geneMod$strand


# should be NULL
sampleFx(0)
sampleFx(100000000)
sampleFx(preMod[lG] + 1)
sampleFx(preMod[1] - 1)


# Done! Time to try more possible ranges...


# Another misc. test
geneData <- list()
geneData$Bob <- list(strand='+', start=200, end=300, exonStart=c(201,221,241), exonEnd=c(210,230,250))

sampleFxBob <- makeMap(geneData$Bob)

sampleFxBob(200)
sampleFxBob(201)
sampleFxBob(202)
sampleFxBob(203)
sampleFxBob(221)
sampleFxBob(249)
sampleFxBob(250)

geneData$Steve <- list(strand='-', start=200, end=300, exonStart=c(201,221,241), exonEnd=c(210,230,250))

sampleFxSte <- makeMap(geneData$Steve)

sampleFxSte(200)
sampleFxSte(201)
sampleFxSte(202)
sampleFxSte(203)
sampleFxSte(221)
sampleFxSte(249)
sampleFxSte(250)




#===================Overall function factory prototype==========================

mapCoord <- function(geneData, dropExonBoundaries = TRUE, silent = FALSE, writeLog = TRUE){
    for (name in names(geneData)) {
        geneData[[name]]$map <- makeMap(geneData[[name]])
    }
    saveRDS(geneData, file = "genedata.rds")

}


# [END]
