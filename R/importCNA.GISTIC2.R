# importCNA.GISTIC2.R
#
#
# Imports copy number aberration data from Firehose database,
# produces an rCNA object saved as a RDS file.
#
#
#'@param fFHG A vector of local fully qualified paths of Firehose GISTIC2 CNA file(s).
#'@param dCNA The local fully qualified path of a directory to store the rCNA RDS file(s).
#'@param silent Boolean. Default: FALSE. Whether or not to write progress information to console.
#'@param noLog Boolean. Default: FALSE. Whether or not to log results.
#
# Opens all files in a requested directory,
# reads them into a data frame, and save the resulting object in RDS format.
#
#'@family importM.COSMIC, importM.TCGA, importM.GISTIC2, importM.ProjectGenie
#
#'@example
# fFHG <- "inst/extdata/devCNA.txt"
# dCNA <- "inst/extdata/dCNA"
# silent <- FALSE
# noLog <- FALSE
#
# importCNA.GISTIC2(fFHG, dCNA, silent = FALSE, noLog = FALSE)
#
# Input - txt (table)
# first 3 columns are gene, NIH locus ID and cytoband, remaining columns are sample numbers
# each row for one gene
# copy number values are in table
#
# Output
# Rownames: HGNC gene symbols.
# Column names: Tumor Sample Barcodes
# Values: CNA as copynumber -2. Floating point number or NA
#
# TODO
# MAKE SURE PATH GIVEN IS VALID AND FILE IS NOT EMPTY/CORRUPT
# READ
# GET RID OF INFORMATION NOT NEEDED (FILTER)
# SAVE AS DATAFRAME AS RDS,  AND STORE IN DIR
# PRINT STATUS
#
#
importCNA.GISTIC2<- function(fFHG, dCNA, silent = FALSE, noLog = FALSE){

    # start log file
    if (!noLog){
        write("Log begin",file="test.log",append=TRUE)
        write(format(Sys.time(), "%a %b %d %X %Y"),file="test.log",append=TRUE)
        write(fFHG,file="test.log",append=TRUE)
        write(dCNA,file="test.log",append=TRUE)
    }

    # check if file exists

    numFiles <- length(fFHG)
    for (file in 1:numFiles){
        if(!file.exists(fFHG[file])){
            sprintf(" ERROR: %s does not exist. Recheck path in fFHG and try again.",fFHG[file] )
            if (!noLog){
                write(fFHG[file],file="test.log",append=TRUE)
                write("ERROR: File did not exist",file="test.log",append=TRUE)
            }
            stop()
        }
    }

    # Update Log

    if (!noLog){
        write("Files exist. Now parsing.",file="test.log",append=TRUE)
    }

    # read and write files in vector


    for (f in 1:numFiles){
        #read file into tab delimited table
        temp <- read.delim(fFHG[f])
        # write to the proper place
        fi <- paste0(dCNA,"/rCNA",f,".rds")
        # get rid of unwanted columns
        #fi$Locus.ID <- NULL
        #fi$Cytoband <- NULL
        # save dataframe
        file.create(fi)
        saveRDS(temp, file=fi)
    }

    # write to console

    if (!silent){
        sprintf("All files in fFHG written into %s. Success! Import from FIREHOSE complete", dCNA)
    }

    # finish log file

    if (!noLog){
        write(format(Sys.time(), "%a %b %d %X %Y"),file="test.log",append=TRUE)
        write("Log end",file="test.log",append=TRUE)
    }


    }

