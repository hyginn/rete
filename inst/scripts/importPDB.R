# importPDB.R
# Non-exported function

library("bio3d")
library("Biostrings")
install.packages("seqnir")

readRDS(MT.RDS)


for (genes in gene_list){
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
        #Perform pairwisealignment..
        ali <- pairwiseAlignment(subjectAA, targetAA)
        # ^^^ Q: what should the gapExtension and gapOpening penalties be?
        for (i in 1:nmatch(ali)){   #For all matched residues:
            i <- index of residue in targetAA
            j <- index of residue in subjectAA
            #^^^ Q: How do I get the indices of residues aligned in target and subject?
            rowNames <- append(rowNames, i)     # rowname will be the index of the residue in target gene
            resIDs <- append(unique(subject$atom[subject$atom$resno==j,]$resid)) #resids are the 3 letter code of the aa
            aaTypes <- append(aaTypes, type of targetAA[i])
            # ^^^ Q: what is aaType? (hydrophobic vs. hydrophilic? or acidic vs. basic? or something else?)

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
        DF <- data.frame("Structure ResID" <- resIDs,
                         "AA Type" <- aaTypes,
                         "Sidechain Centroid" <- sidechainXYZ,
                         "C.alpha Centroid" <- calphaXYZ)
        rownames(DF) <- rowNames

        # append it to the results list
        results <- append(results, DF)
    }
    #store...
    store results

}
