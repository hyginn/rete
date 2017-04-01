# importPDB.R
# Non-exported function

#=============================STOP & SAVE======================================
# Provide a mechanism to stop:
# REFERENCE:
#   http://stackoverflow.com/questions/38755283/break-loop-with-keyboard-input-r
# requires library(tcltk2)
win1 <- tktoplevel()
butStop <- tkbutton(win1, text = "Stop",
                    command = function() {
                        assign("stoploop", TRUE, envir = .GlobalEnv)
                        tkdestroy(win1)
                    })

butSave <- tkbutton(win1, text = "Pause & Save",
                    command = function() {
                        assign("stoploop", TRUE, envir = .GlobalEnv)
                        print("Finishing the current gene..")
                        saveRDS(genes, file = "inst/tmp/genes.rds")
                        saveRDS(GeneData, "inst/extdata/GeneData.rds")
                        tkdestroy(win1)
                    })
tkgrid(butStop)
tkgrid(butSave)
stoploop <- FALSE

readRDS(MT.rds)

#=============================(LOAD & CONTINUE)/START==========================
# Check if there was a saved progress from before.
if (file.exists("inst/tmp/genes.rds")){
    genes <- readRDS("inst/tmp/genes.rds")

} else {
    genes <- list(gene names from MT)
}

while (length(genes > 0) && !stoploop){
    results <- list()

    #=========================GET TARGET SEQUENCE==============================

    # The current gene is always the first element of the list.
    # This element will be removed from the list of genes at the end of each
    # iteration of the while loop.This is necessary for pausing the process and
    # continuing it later.

    curr_gene <- genes[[1]]

    # Translate the cDNA of the gene to its aa sequence
    targetAA <- AAString(paste(a.a. sequence of the gene, collapse = "")

    for (pdbID_chainID matched to gene){
        pdbID <- pdbid retrieved from the MT object
        chainID <- chainid retrieved from the MT object

        #=====================READ SUBJECT PDB FILE============================

        # download the PDB file.
        subject <- read.pdb(pdbID)
        # extract the desired chain and get its aa sequence.
        chainAAIndices <- atom.select(subject,
                                      "calpha",
                                      chain = chainID,
                                      value = TRUE)
        subjectAA <- AAString(paste(aa321(chainAAIndices$atom$resid), collapse = ""))

        #=====================ALIGNMENT========================================

        # retrieve the start and end indices for the target & subject sequences:
        start_target <- the starting index of the alignment in targetAA from MT
        end_target <- the last index of the alignment in targetAA from MT
        start_subject <- the starting index of the alignment subjectAA from MT
        end_subject <- the last index of the alignment in target AA from MT

        # perform pairwise alignment with default parameters.
        ali <- pairwiseAlignment(subjectAA[start_subject, end_subject],
                                 targetAA[start_target, end_target])

        #======================================================================

        # create a data frame with nmatch(ali) rows.
        DF <- data.frame("ResID" = c(1:nmatch(ali)),
                         "AA Type" = c(1:nmatch(ali)),
                         "Sidechain Centroid" = c(1:nmatch(ali)),
                         "C.alpha Centroid" = c(1:nmatch(ali)))

        # for all matched residues:
        for (n in 1:nmatch(ali)){
            i <- index of residue in targetAA
            j <- index of residue in subjectAA

            # rowname will be the index of the residue in target gene:
            rownames(DF)[n] <- i

            #=================RESIDUE ID & AA TYPE=============================

            # resids are the 3 letter code of the aa:
            resID <- unique(subject$atom[subject$atom$resno==j,]$resid)
            # type of AA..
            aaType <- type of targetAA[i]
            # ^^^ Q: what does it mean by amino acid type?
            #        (hydrophobic vs. hydrophilic? or acidic vs. basic? or
            #        something else?) Still need to figure out why we need
            #        these.. They are mentioned in the design task page of
            #        the GeneData object.

            #=================SIDECHAIN CENTROID===============================

            # subset the sidechain atoms; do not include hydrogen atoms:
            sidechainInd <- combine.select(atom.select(subject,
                                                       "sidechain",
                                                       chain = chainID,
                                                       resno = j),
                                           atom.select(subject,
                                                       "noh")
            sideChain <- subject$atom[sidechainInd$atom,]
            # calculate the sidechain centroid:
            sidechainXYZ <- list(x = mean(sideChain$atom$x),
                                 y = mean(sideChain$atom$y),
                                 z = mean(sideChain$atom$z))

            #=================ALPHA CARBON CENTROID============================

            # subset the alpha carbon:
            C.alpha <- atom.select(subject,
                                   "calpha",
                                   chain = chainID,
                                   resno = j)
            # get the x, y, z coordinates of the alpha carbon:
            C.alphaXYZ <- list(x = C.alpha$atom$x,
                               y = C.alpha$atom$y,
                               z = C.alpha$atom$z)

            #==================================================================

            # append the data to the data frame
            DF[n,] <- c(resID, aaType, sidechainXYZ, C.alphaXYZ)
        }

        # add the data frame to the list of data frames in GeneData$3D.
        # the index in the list will be of format <pdbID_chainID>.
        GeneData$3D[paste0(PDB_ID, "_", chain_ID)] <- DF
    }

    # remove the current gene from the list of genes.
    genes[[1]] <- NULL
}

# save the GeneData object:
saveRDS(GeneData, "inst/extdata/GeneData.rds")

# Remove the genes.rds file (if it exists) after 3D coordinates of
# all genes have been populated.
if (file.exists("inst/tmp/genes.rds")){
    file.remove("inst/tmp/genes.rds")
}
