# testCombine.R

context("test for Combine SNV and CNA files")

# ==== BEGIN SETUP AND PREPARE =================================================
# set up tempdir and tempfiles for output files and to be deleted after tests
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
# ==== END SETUP AND PREPARE ===================================================
# setting up small input rds data
# building SNV.rds with 2 genes
testSNV = tempfile(fileext="rds")
tmpSNV <- as.data.frame(matrix(nrow=2, ncol=12))
rownames(tmpSNV) <- c("PTCHD2", "NR113")
colnames(tmpSNV) <- c("sym", "chr",	"start", "end",	"strand", "class", "type", "aRef", "a1", "a2", "TCAG", "UUID")
# only entering data useful for combine function. e.g. "sym". "start", "class"
tmpSNV["PTCHD2", "sym"] <- "PTCHD2"
tmpSNV["PTCHD2", "start"] <- 11534434
tmpSNV["PTCHD2", "class"] <- "Missense_Mutation"
tmpSNV["NR113", "sym"] <- "NR113"
tmpSNV["NR113", "start"] <- 161230879
tmpSNV["NR113", "class"] <- "Missense_Mutation"
saveRDS(tmpSNV, file=testSNV)

#building CNA.rds with 2 genes with 2 tumor type
testCNA = tempfile(fileext="rds")
tmpCNA <- as.data.frame(matrix(nrow=2, ncol=2))
rownames(tmpCNA) <- c("ACAP3", "ACTRT2")
colnames(tmpCNA) <- c("TCGA-OR-A5J1-01A-11D-A29H-01", "TCGA-OR-A5J2-01A-11D-A29H-01")
tmpCNA["ACAP3", "TCGA-OR-A5J1-01A-11D-A29H-01"] <- 0.030
tmpCNA["ACAP3", "TCGA-OR-A5J2-01A-11D-A29H-01"] <- -0.070
tmpCNA["ACTRT2", "TCGA-OR-A5J1-01A-11D-A29H-01"] <- 0.030
tmpCNA["ACTRT2", "TCGA-OR-A5J2-01A-11D-A29H-01"] <- 0.070
saveRDS(tmpCNA, file=testCNA)

test_that("parameter errors are correctly handled", {
    # try no parameter input
    expect_error(combineSNV_CNA(writeLog=FALSE),
                 'Input vectors are empty')
    # sane input vector
    input_vector = c(testSNV, testCNA)
    # try non existing input file for SNV input vector
    expect_error(combineSNV_CNA(c("nonExistent.rds"), input_vector[2], writeLog=FALSE),
                 '.checkArgs> "fSNV" error: file "nonExistent.rds" does not exist')

    # try non existing input file for CNA input vector
    expect_error(combineSNV_CNA(input_vector[1], c("nonExistent.rds"), writeLog=FALSE),
                 '.checkArgs> "fCNA" error: file "nonExistent.rds" does not exist')

    # try invalid output fgX file name (e.g. wrong extension)
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], "out.txt", writeLog=FALSE),
                 'Output file has to be .rds extension!')

    # try null output fgX file name
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], NULL, writeLog=FALSE),
                 'argument is of length zero')

    # try null silent option
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], silent=NULL, writeLog=FALSE),
                 '.checkArgs> "silent" mode error')

    # try null writeLog option
    expect_error(combineSNV_CNA(input_vector[1], input_vector[2], writeLog=NULL),
                 '.checkArgs> "writeLog" mode error')

    # check if log file has been written
    expect_false(file.exists(logName))
})


test_that("a sane input gives an expected output", {
    # provide one valid SNV and one valid CNA file with small number of genes as input

    input_vector = c(testSNV, testCNA)
    testF <- tempfile(fileext = ".rds")
    # test if a log file is saved if log is True
    combineSNV_CNA(input_vector[1], input_vector[2], fgX=testF, silent=TRUE, writeLog=TRUE)
    expect_true(file.exists(logName))
    # correct log file
    thisLog <- readLines(logName)
    expect_true(grepl("event \\| title  \\| combineSNV_CNA",       thisLog[1]))
    expect_true(grepl("event \\| time",                            thisLog[2]))

        # call
    expect_true(grepl("^event \\| call   \\| combineSNV_CNA\\(",   thisLog[3]))
    expect_true(grepl(paste0("fNameSNV = \"", input_vector[1], "\""),    thisLog[3]))
    expect_true(grepl(paste0("fNameCNA = \"", input_vector[2], "\""),    thisLog[3]))
    expect_true(grepl(sprintf("fgX = \"%s\"", testF),                  thisLog[3]))
    expect_true(grepl("silent = TRUE",                             thisLog[3]))
    expect_true(grepl("writeLog = TRUE",                           thisLog[3]))
    expect_true(grepl(")$",                                        thisLog[3]))

    expect_true(grepl("event | note   | Combine 10 samples in total.",  thisLog[4]))
    expect_true(grepl("event | end",                               thisLog[5]))
    expect_true(grepl("^$",                                        thisLog[6]))


    # test if output file fgX is saved under the path
    expect_true(file.exists(testF))
    # remove tempfile
    expect_true(file.remove(testF))
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
# remove testCNA and testSNV
file.remove(testCNA)
file.remove(testSNV)
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
