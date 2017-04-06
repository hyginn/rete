# importPDB.R

#' Populate a GeneData object with the 3D coordinates of each gene.
#'
#' \code{importPDB} Recieves a data frame with rownames being PDB-ID/chain, and
#'                  columns being the first and last residue for which 3D
#'                  coordinates have been mapped.
#'
#' @section Read pattern PDB file:
#'   Reads a PDB file using bio3d::read.pdb(), retrieves the actual sequence
#'   from the coordinate record (column $atom), and converts this sequence
#'   into an AAstring of 1-letter amino acid codes.
#'
#'  @section Obtain start & end indices:
#'    Retrieves the start and end indices of the alignment for the subject
#'    abd pattern sequences from GeneData$xyzMap
#'
#' @section Populate GeneData object:
#'   Fetches data from PDB files and creates a data frame to be attached to
#'   GeneData$xyz, with one row for each residue matched. This section includes
#'   the following subsections:
#'
#'   @section Residue ID & AA type:
#'     Fetches the residue ID and amino acid type of the each residue matched.
#'     residue ID is the position of the amino acid within the actual sequence
#'     and amino acid type is the 1-letter amino acid code of the residue.
#'
#'   @section Sidechain centroid:
#'     Calculates the sidechain centroid of each residue by fetching and
#'     averaging over all x,y,z coordinates of atoms making up the sidechain of
#'     the residue matched in the pattern sequence.
#'
#'   @section Alpha carbon centroid:
#'     Calculates the alpha carbon centroid of each residue by fecthing the
#'     x, y, z coordinates of the alpha carbon of the residue matched in the
#'     pattern sequence.
#'
#' @param XYZMap the $xyzMap column of a GeneData object that contains PDB-ID/chain of
#'                  aligned sequences and the the first and last residue for which 3D
#'                  coordinates have been mapped.
#' @param silent logical. Whether output will be to console will be suppressed.
#'                  Default FALSE.
#' @param writeLog logical. Whether an event entry for the log file is to be
#'                  written. Default TRUE.
#' @return <description>.
#'
#' @family <optional description of family>
#'
#' @return N/A. This function is invoked for its side effect of writing 3D
#'   coordinates to a GeneData object.
#'
#' @seealso \code{\link{SEL3D}} Blasts input sequence and return MT object with
#'   match information.
#'
#' @examples
#' importPDB(myGeneData)

importPDB <- function(XYZMap, silent = FALSE, writeLog = TRUE){
    # Check input sanity:
    if (mode(XYZMap) != "list" ||
        mode(silent) != "logical" ||
        mode(writeLog) != "logical"){
        stop("Invalid input!")
    }

    if (!silent){ cat("Populating 3D coordinates of a gene...\n") }
    # For pdb.IDs matched to this gene:
    for (row in rownames(XYZMap)){
        if (!silent){ cat("\tFetching data for", XYZMap[row, "pdb.id"], "\n") }
        # Get the pdb.id and chain.id of this match.
        pdbIDchainID <- strsplit(XYZMap[row, "pdb.id"], "_")[[1]]
        pdbID <- pdbIDchainID[1]
        chainID <- pdbIDchainID[2]

        #=====================READ PATTERN PDB FILE============================
        # Download the PDB file.
        pattern <- bio3d::read.pdb(pdbID, verbose = FALSE)
        # Extract the desired chain and get its aa sequence.
        chainAAIndices <- atom.select(pattern,
                                      "calpha",
                                      chain = chainID,
                                      value = TRUE)
        patternAA <- AAString(paste(aa321(chainAAIndices$atom$resid),
                                            collapse = ""))


        #**********************************************************************
        #* From this point on, the code is not tested for its accuracy or     *
        #* results and needs further improvements/testing/validation.         *
        #**********************************************************************

        #=====================OBTAIN START & END INDICES=======================
        # Retrieve the start (s) and end (e) indices for the pattern &
        #   subject sequences:
        pStart <- XYZMap$q.start
        pEnd <- XYZMap$q.end
        sStart <- XYZMap$s.start
        sEnd <- XYZMap$s.end

        #=====================POPULATE GENEDATA OBJECT=========================
        # Create a data frame with number of matched residues rows.
        nMatch <- pEnd - pStart + 1
        DF <- data.frame("ResID" = numeric(nMatch),
                         "AA Type" = numeric(nMatch),
                         "Sidechain Centroid" = list(x = numeric(nMatch),
                                                     y = numeric(nMatch),
                                                     z = numeric(nMatch)),
                         "C.alpha Centroid" = list(x = numeric(nMatch),
                                                   y = numeric(nMatch),
                                                   z = numeric(nMatch)),
                         stringsAsFactors = FALSE)

        # For all matched residues:
        for (n in 1:nMatch){
            if (!silent){
                cat("\t\tFetching data for macthing residue", n, "\n")
            }
            # Find the position (residue number) of the aligned residues in
            # the pattern & subject sequences.
            pResNo <- pStart + n - 1
            sResNo <- sStart + n - 1

            # Row name will be the index of the residue in the subject:
            #   - I had to add "res" to the beginning of each rowname since
            #     R was giving me an error for having duplicate row names:
            #    -> Error in `row.names<-.data.frame`(`*tmp*`, value = value):
            #       duplicate 'row.names' are not allowed"
            rownames(DF)[n] <- paste0("res", sResNo)

            #=================RESIDUE ID & AA TYPE=============================
            # Resid is the index of the residue in the subject:
            resID <- sResNo
            # Amino acid types are the 1-letter code of the aa:
            aaType <- toString(patternAA[pResNo])

            #=================SIDECHAIN CENTROID===============================
            # First, we need to find the {residue number + insert code} of the
            #   residue in the PDB file. To do this, we create a vector of all
            #   _resno's_ and _insertion codes_ in the PDB file. Then we simply
            #   subset the pResNo value from the two vector.
            #     - This precedure assumes that there are no missing PDB data
            #       for any residue.
            #     - {value = True} causes _atom.select_ to return a pdb object
            #       instead of indices.
            pdbResNoVector <- atom.select(pattern,
                                          "calpha",
                                          value=TRUE,
                                          verbose = FALSE)$atom$resno
            pdbinsertVector <- atom.select(pattern,
                                           "calpha",
                                           value=TRUE,
                                           verbose = FALSE)$atom$insert
            pdbResNo <- pdbResNoVector[pResNo]
            pdbInsert <- pdbinsertVector[pResNo]

            # Subset the sidechain atoms; do not include hydrogen atoms:
            #   - sidechainInd$atom is a vector of indices.
            sidechainInd <- bio3d::combine.select(bio3d::atom.select(pattern,
                                                                     "sidechain",
                                                                     chain = chainID,
                                                                     resno = pdbResNo,
                                                                     insert = pdbInsert),
                                                  bio3d::atom.select(pattern,
                                                                     "noh"),
                                                  verbose = FALSE)
            sidechain <- pattern$atom[sidechainInd$atom,]
            # calculate the sidechain centroid:
            sidechainXYZ <- list(x = mean(sideChain$x),
                                 y = mean(sideChain$y),
                                 z = mean(sideChain$z))

            #=================ALPHA CARBON CENTROID============================
            # subset the alpha carbon:
            #  - {value = True} causes _atom.select_ to return a pdb object
            #    instead of indices.
            C.alpha <- bio3d::atom.select(pattern,
                                          "calpha",
                                          chain = chainID,
                                          resno = pdbResNo,
                                          insert = pdbInsert,
                                          value = TRUE)
            # get the x, y, z coordinates of the alpha carbon:
            C.alphaXYZ <- list(x = C.alpha$atom$x,
                               y = C.alpha$atom$y,
                               z = C.alpha$atom$z)

            #==================================================================
            # Place the data in the data frame
            DF[n,] <- c(resID, aaType, sidechainXYZ, C.alphaXYZ)
        }

        # Add the data frame to the list of data frames in GeneData$3D.
        # the index in the list will be of format <pdbID_chainID>.
        GeneData$xyz[paste0(PDB_ID, "_", chain_ID)] <- DF
    }

    if(writeLog) {
        myTitle <- "importPDB"

        # Compile function call record
        myCall <- character()
        myCall[1] <- "importPDB("
        myCall[2] <- sprintf("XYZMap = \"%s\", ", XYZMap)
        myCall[3] <- sprintf("silent = %s, ", as.character(silent))
        myCall[4] <- sprintf("writeLog = %s)", as.character(writeLog))
        myCall <- paste0(myCall, collapse = "")

        # Record progress information
        myNotes <- character()
        myNotes <- c(myNotes, sprintf("Read xyzMap from file: %s", XYZMap))
        myNotes <- c(myNotes, sprintf("Populated 3D coordinates for %s genes",
                                      length(GeneData$xyz)))

        # indicate output object name(s)
        myOutput = c("GeneData")

        # send info to log file
        logEvent(eventTitle = myTitle,
                 eventCall = myCall,
                 notes = myNotes,
                 output = myOutput)
    }
}

#=============================SCRIPT===========================================
#=============================STOP & SAVE======================================
# Provide a mechanism to stop:
# REFERENCE:
#   http://stackoverflow.com/questions/38755283/break-loop-with-keyboard-input-r
# Requires library(tcltk2)
PROGRESSFILE <- "inst/tmp/PROGRESSFILE.rds"
win1 <- tktoplevel()
butStop <- tkbutton(win1, text = "Stop",
                    command = function() {
                        assign("stoploop", TRUE, envir = .GlobalEnv)
                        tkdestroy(win1)
                    })

butSave <- tkbutton(win1, text = "Pause & Save",
                    command = function() {
                        assign("stoploop", TRUE, envir = .GlobalEnv)
                        print("Finishing the current gene...")
                        saveRDS(genes, file = PROGRESSFILE)
                        print("A progress file has been created under inst/tmp/")
                        saveRDS(GeneData, "inst/extdata/GeneData.rds")
                        tkdestroy(win1)
                    })

tkgrid(butStop)
tkgrid(butSave)
stoploop <- FALSE

# Load the GeneData object.
myGeneData <- readRDS("inst/extdata/GeneData.rds")

#=============================(LOAD & CONTINUE)/START==========================
# Check if there exists a saved progress file from before.
if (file.exists("inst/tmp/genes.rds")){
    genes <- read.RDS(PROGRESSFILE)
} else {
    genes <- names(myGeneData)
}

#=============================RUN IMPORTPDB()==================================
while (length(genes > 0) && !stoploop){
    # The current gene is always the first element of the list.
    # This element will be removed from the list of genes at the end of
    #   each iteration of the while loop. This is necessary for pausing
    #   the process and continuing it later.
    currentGene <- genes[[1]]
    XYZMap <- myGeneData[[currentGene]]$xyzMap

    importPDB(XYZMap)

    # Remove the current gene from the list of genes.
    genes[[1]] <- NULL
}

# Save the GeneData object:
saveRDS(GeneData, "inst/extdata/GeneData.rds")

# Remove the genes.rds file (if it exists) after 3D coordinates of
# all genes have been populated.
if (file.exists(PROGRESSFILE)){ file.remove(PROGRESSFILE) }


