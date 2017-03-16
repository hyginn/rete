# importCNA.GISTIC2.R
#
#
#'@function Imports copy number aberration data from Firehose database,
# produces an rCNA object saved as a RDS file.
#
#
#'@param fFHG Local fully qualified path of Firehose GISTIC2 CNA file.
#'@param dCNA The local fully qualified path of a directory to store the rCNA RDS file(s).
#'@param silent Boolean. Default: FALSE. Whether or not to write progress information to console.
#'@param noLog Boolean. Default: FALSE. Whether or not to log results.
#
#
#'@family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#
#'@example
# fFHG <- "inst/extdata/devCNA.txt"
# dCNA <- "inst/extdata/dCNA"
# silent <- FALSE
# writeLog <- TRUE
#
# importCNA.GISTIC2(fFHG, dCNA, silent = FALSE, writeLog = TRUE)
#
# Input - all_data_by_genes.txt from the GISTIC2 output from FireHose
# saved locally
#
# Output - rCNA file
#
#
#

importCNA.GISTIC2<- function(fFHG, dCNA, silent = FALSE, writeLog = TRUE){

    # check if file exists

    if(!file.exists(fFHG)){
            stop(fFHG)
    }

    # parameter check

    cR <- character()
    cR <- c(cR, .checkArgs(fFHG,        like = "devCNA.txt",   checkSize = TRUE))
    cR <- c(cR, .checkArgs(dCNA,        like = "inst/extdata/dCNA",
                           checkSize = TRUE))
    cR <- c(cR, .checkArgs(silent,       like = logical(1), checkSize = TRUE))
    cR <- c(cR, .checkArgs(writeLog,     like = logical(1), checkSize = TRUE))

    if(length(cR) > 0) {
        stop(cR)
    }

    # read and write files in vector
    # commented out part that would handle multiple files at one time

    #for (f in 1:numFiles){
        #read file into tab delimited table
        #temp <- read.delim(fFHG)
        temp <- readr::read_delim(file = fFHG,
                                  delim = "\t", col_names =TRUE,
                                  col_types = c(col_character(),col_integer(),
                                                col_character(), col_double()))
        # check file
        if (colnames(temp)[1] != "Gene Symbol"){
            stop("Not a valid file")
        }
        if (colnames(temp)[2] != "Locus ID"){
            stop("Not a valid file")
        }
        if (colnames(temp)[3] != "Cytoband"){
            stop("Not a valid file")
        }

        # write to the proper place so it doesn't get overwritten
        namefile<- basename(fFHG)
        fi <- paste0(dCNA, "/", namefile, ".rds")

        # get rid of unwanted columns
        #fi$Locus.ID <- NULL
        #fi$Cytoband <- NULL

        # save dataframe
        file.create(fi)
        saveRDS(temp, file=fi)
    #}

    # write to console

    if (silent){
        sprintf("File written into %s. Success! Import from FIREHOSE complete", dCNA)
    }


# ==== WRITE LOG ===========================================================
    if(writeLog) {
        msg    <- sprintf("event > %s from importNet.STRING()\n", Sys.time())
        msg[2] <- "note >   Returned gG object with"
        msg[3] <- sprintf("%d vertices and %d edges.\n",
                              igraph::gorder(gG),
                              igraph::gsize(gG))
        msg[4] <- sprintf("out > %s", "<object UUID and attributes>")
        msg[5] <- "\n"
        logMessage(msg)
    }

        return(gG)

    }

