# importPDB.R
# Non-exported function

library("bio3d")
library("Biostrings")
install.packages("seqnir")

readRDS(MT.RDS)

# Provide a mechanism to stop:
# REFERENCE:
#   http://stackoverflow.com/questions/38755283/break-loop-with-keyboard-input-r
library(tcltk)
win1 <- tktoplevel()
butStop <- tkbutton(win1, text = "Stop",
                    command = function() {
                        assign("stoploop", TRUE, envir=.GlobalEnv)
                        tkdestroy(win1)
                    })

butSave <- tkbutton(win1, text = "Pause & Save",
                    command = function() {
                        assign("stoploop", TRUE, envir=.GlobalEnv)
                        print("Finishing the current gene..")
                        saveRDS(genes, file = "inst/tmp/genes.rds")
                        tkdestroy(win1)
                    })
tkgrid(butStop)
tkgrid(butSave)
stoploop <- FALSE

if (file.exists("inst/tmp/genes.rds")){
    genes<- readRDS("inst/tmp/genes.rds")

} else {
    genes <- list(gene symbols from MT)
}

while (length(genes > 0) && !stoploop){
    # The current gene is always the first element of the list.
    # This element will be removed from the list of genes at the end of each
    # iteration of the while loop.This is necessary for pausing the process and
    # continuing it later.
    curr_gene <- genes[[1]]

    results <- list()
    # Translate the gene to its aa sequence
    targetAA <- AAString(paste(a.a. sequence of the gene, collapse="")

    for (PDB-ID/chains matched to gene){
        # vectors for the data.frame that will be created at the end of this for loop.
        rowNames <- c()
        resIDs <- c()
        aaType <- c()
        sidechainXYZ <- c()
        calphaXYZ <- c()

        # load the PDB file.
        subject <- read.pdb(PDB-ID)
        # get the aa sequence of the subject
        subjectAA <- AAString(paste(PDB$seqres, collapse=""))
        # retrieve the start and end indices for the target and subject sequences.
        start_target <- the starting index of the alignment in targetAA from MT
        end_target <- the last index of the alignment in targetAA from MT
        start_subject <- the starting index of the alignment subjectAA from MT
        end_subject <- the last index of the alignment in target AA from MT
        # Perform pairwisealignment..
        ali <- pairwiseAlignment(subjectAA[start_subject, end_subject],
                                 targetAA[start_target, end_target])
        # ^^^ Q: what should the gapExtension and gapOpening penalties be?

        for (i in 1:nmatch(ali)){   #For all matched residues:
            i <- index of residue in targetAA
            j <- index of residue in subjectAA
            #^^^ Q: How do I get the indices of residues aligned in target and subject?
            rowNames <- append(rowNames, i)     # rowname will be the index of the residue in target gene
            resIDs <- append(unique(subject$atom[subject$atom$resno==j,]$resid)) #resids are the 3 letter code of the aa
            aaTypes <- append(aaTypes, type of targetAA[i])
            # ^^^ Q: what does it mean by amino acid type?
            #        (hydrophobic vs. hydrophilic? or acidic vs. basic? or something else?)

            # specify which atoms in subject$atom$elety are for alpha carbon, COOH and NH2.
            # The rest of the atoms must be in the side chain
            nonSideChain_elety <- c("N", "C", "O", "H", "CA")
            thisResData <- subject$atom[subject$atom$resno== j,]

            #If there is data for this residue:
            if (nrow(thisResData)!=0){
                #x, y, z coordinates of the side chain atoms
                sideChainXValues <- thisResData[!(thisResData$elety %in% nonSideChain_elety), "x"]
                sideChainYValues <- thisResData[!(thisResData$elety %in% nonSideChain_elety), "y"]
                sideChainZValues <- thisResData[!(thisResData$elety %in% nonSideChain_elety), "z"]

                # calculate the average of the x, y, z coordinates
                avgX <- sum(sideChainXValues)/length(sideChainXValues)
                avgY <- sum(sideChainYValues)/length(sideChainYValues)
                avgZ <- sum(sideChainZValues)/length(sideChainZValues)
                sidechainXYZ <- append(sidechainXYZ, c(avgX, avgY, avgZ))

                # x, y, z coordinates of the alpha carbon
                calphaXYZ <- append(calphaXYZ, c(thisResData[thisResData$elety=="CA", "x"),
                                                 thisResData[thisResData$elety=="CA", "y"],
                                                 thisResData[thisResData$elety=="CA", "z"])
            }
        }

        # create the data frame.
        DF <- data.frame("ResID" <- resIDs,
                         "AA Type" <- aaTypes,
                         "Sidechain Centroid" <- sidechainXYZ,
                         "C.alpha Centroid" <- calphaXYZ)
        rownames(DF) <- rowNames

        # append it to the results list
        results <- append(results, DF)
    }
    #store results...
    GeneData$3D[paste0(PDB_ID, chain_ID)] <- results

    # Remove the current gene from the list of genes.
    genes[[1]] <- NULL
}

# Remove the genes.rds file (if it exists) after 3D coordinates of
# all genes have been populated.
if (file.exists("inst/tmp/genes.rds")){
    file.remove("inst/tmp/genes.rds")
}
