# importCNA.GISTIC2.R
#
#
#'@function Imports copy number aberration data from Firehose database,
# produces an rCNA object saved as a RDS file.
#
# ########## ToDo (bs): Fix Roxygen comments
# ########## ToDo (bs): Add Details section
#
#'@param fFHG Local fully qualified path of Firehose GISTIC2 CNA file.
#'@param dCNA The local fully qualified path of a directory to store the rCNA RDS file(s).
#'@param silent Boolean. Default: FALSE. Whether or not to write progress information to console.
#'@param noLog Boolean. Default: FALSE. Whether or not to log results.
########## ToDo (bs):  replace noLog with writeLog
#'@return gG object
########## ToDo (bs):  not a return value
#'@notRun {
########## ToDo (bs):  Not run needs to be part of the Examples section
########## ToDo (bs):  Why is family in the notRun block? 
#'@family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#}
#'@example
# fFHG <- "inst/extdata/devCNA.txt"
# dCNA <- "inst/extdata/dCNA"
# silent <- FALSE
# writeLog <- TRUE
# importCNA.GISTIC2("inst/extdata/devCNA.txt", "inst/extdata/dCNA")
#'@export
#
#
#

importCNA.GISTIC2<- function(fFHG, dCNA, silent = FALSE, writeLog = TRUE){

    # check if file exists

    if(!file.exists(fFHG)){
########## ToDo (bs):  Should be handled by .checkArgs()
            stop(fFHG)
    }

    # parameter check

    cR <- character()
    cR <- c(cR, .checkArgs(fFHG,        like = "devCNA.txt",   checkSize = TRUE))
########## ToDo (bs):  wrong check: need to check for FILE_E
    cR <- c(cR, .checkArgs(dCNA, like = "DIR", checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # read and write files in vector
    # commented out part that would handle multiple files at one time
########## ToDo (bs):  Why? The spec. mentions multiple files.
########## ToDo (bs):  But I agree that this is probably not a great idea.
########## ToDo (bs):  Still, this comment should be removed.


    temp <- readr::read_delim(file = fFHG,
                        delim = "\t", col_names =TRUE,
                        col_types = c(col_character(),col_integer(),
                        col_character(), col_double()))
########## ToDo (bs):  Comment: what are these columns?
        # check file
        if (colnames(temp)[1] != "Gene Symbol"){
            stop("Not a valid file")
########## ToDo (bs):  This is not an informative error message
        }
        if (colnames(temp)[2] != "Locus ID"){
########## ToDo (bs):  Is this guarenteed to be an integer? Do we need this is all?
            stop("Not a valid file")
        }
        if (colnames(temp)[3] != "Cytoband"){
########## ToDo (bs):  Do we need this column?
            stop("Not a valid file")
        }

        # write to the proper place so it doesn't get overwritten
        namefile<- basename(fFHG)
########## ToDo (bs):  Misleading variable name. Use CamelCase.
        fi <- paste0(dCNA, "/", namefile, ".rds")
########## ToDo (bs):  Wrong. Need to use file.path() for platform independence.

        # get rid of unwanted columns
########## ToDo (bs):  Right - but these should be skipped while reading.
        temp$`Locus ID` <- NULL
        temp$Cytoband <- NULL

        # save dataframe
        file.create(fi)
########## ToDo (bs):  Unnecessary
        saveRDS(temp, file=fi)
    #}
########## ToDo (bs):  Why is there a commented-out brace in your submitted code?

    # write to console

    if (silent){
########## ToDo (bs):  That has to be _not_ silent!
########## ToDo (bs):  This is a bug. It should have been caught while testing.
########## ToDo (bs):  Need to update tests.
        sprintf("File written into %s. Success! Import from FIREHOSE complete", dCNA)
########## ToDo (bs):  Need to record how many samples, how many genes.
    }


    # call df2gG

    .df2gG(fi,call = "import CNA from GISTIC2")
########## ToDo (bs):  What. Is. This? 

########## ToDo (bs):  Missing notes
########## ToDo (bs):  Missing metadata    
    
    
# ==== WRITE LOG ===========================================================
    # send info to log file
    logEvent(eventTitle = myTitle,
             eventCall = myCall,
             #                input = character(),
             notes = myNotes,
             output = myOutput)
########## ToDo (bs):  These parameters have not been defined.


        return(gG)
########## ToDo (bs):  Error: No return value!

}
# END
