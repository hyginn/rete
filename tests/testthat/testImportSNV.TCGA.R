# testImportSNV.TCGA.R

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================

#context("Tools for importing MAF files into a rSNV file")
options(stringsAsFactors=FALSE)

test_that("parameter errors are correctly handled", {
    # fMAF is not provided
    expect_error(importSNV.TCGA())
    # fMAF is NULL/empty vector
    expect_error(importSNV.TCGA(fMAF = NULL))
    # fMAF input vector contains non-string values
    expect_error(importSNV.TCGA(fMAF = 123))
    # file in fMAF could not be opened/does not exist
    expect_error(importSNV.TCGA(fMAF = c("nosuchfileMAF.maf")))

    # if rSNV parameter is provided, rSNV file does not exist
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = "nosuchfileSNV"))
    # if rSNV parameter is provided, but is NULL
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = NULL))
    # if rSNV parameter is provided, but is non-string value
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = 123))
    # if rSNV is given vector input
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = c("abc", "xyz")))
    # if rSNV is not provided, new file could not be created (file name conflict)
    expect_error(importSNV.TCGA(fMAF = testMAF))

    # if na.rm is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,na.rm = 1))
    expert_error(importSNV.TCGA(fMAF = testMAF,na.rm = NULL))
    expert_error(importSNV.TCGA(fMAF = testMAF,na.rm = logical()))
    expert_error(importSNV.TCGA(fMAF = testMAF,na.rm = c(TRUE, TRUE)))

    # if silent is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,silent = 1))
    expert_error(importSNV.TCGA(fMAF = testMAF,silent = NULL))
    expert_error(importSNV.TCGA(fMAF = testMAF,silent = logical()))
    expert_error(importSNV.TCGA(fMAF = testMAF,silent = c(TRUE, TRUE)))

    # write.log is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,writeLog = 1))
    expert_error(importSNV.TCGA(fMAF = testMAF,writeLog = NULL))
    expert_error(importSNV.TCGA(fMAF = testMAF,wrtieLog = logical()))
    expert_error(importSNV.TCGA(fMAF = testMAF,writeLog = c(TRUE, TRUE)))
})

test_that("a sane input gives an expected output with differnt input parameters", {

    # make small input (testIn)
    # have an input file, rows too big to import manually here (testIn)

    # make small output for that input manually (testOut)
    write("# metadata \n# 1 \n# 2  \n# 3 \n# 4 \n# 5 \n# metadata END\n", "testOut")
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")
    write(header, "testOut", ncolumns = 12, sep = "\t", append = TRUE)
    sample1 <-c("TMEM201", "chr1", 9607661, 9607661, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("PTCHD2", "chr1", 11516130, 11516130, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA.QQ.A8VF.01A.11D.A37C.09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample4 <- c("TXLNA", "chr1",	32192719,	32192719,	"+",	"Silent",	"SNP",	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample5 <- c("NFYC", "chr1",	40770537,	40770537,	"+",	"Silent",	"SNP",	"C",	"C",	"A", "TCGA-QQ-A8VF-01A-11D-A37C-09",	"2109690d-e8b1-4b5c-9a39-21af273324cf")

    trueOut <- rbind(sample1,sample2,sample3,sample4,sample5)

    # run the module on only fMAF as parameter (moduleOut)
    moduleOut <- importSNV.TCGA(c("testIn"))
    # check if manual output = module output
    # skip 7 for metadata lines just for now
    module <- read.csv("moduleOut", skip=7, sep = "\t")
    all.equal(trueOut, module)

    # test input with rSNV as input
    #  have an input file, rows too big to import manually here (testInPre)

    # test output for rSNV as input
    sample1 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("TXLNA", "chr1",	32192719,	32192719,	"+",	"Silent",	"SNP",	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("NFYC", "chr1",	40770537,	40770537,	"+",	"Silent",	"SNP",	"C",	"C",	"A", "TCGA-QQ-A8VF-01A-11D-A37C-09",	"2109690d-e8b1-4b5c-9a39-21af273324cf")
    trueOutPre <- rbind(sample1,sample2,sample3)

    # run the module on fMAF and rSNV as parameters
    moduleOutPre <- importSNV.TCGA(c("testInPre"), "module02")
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # run the module on fMAF and rSNV and silent as parameters
    moduleOutPre <- importSNV.TCGA(c("testInPre"), "module02", silent=TRUE)
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # run the module on fMAF, rSNV, silent and writeLog as parameters
    moduleOutPre <- importSNV.TCGA(c("testInPre"), "module02", silent=FALSE, writeLog=TRUE)
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # test input with rSNV as input and NA values
    #  have an input file, rows too big to import manually here (testInNA)

    # test output for rSNV as input and NA values
    sample1 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")

    trueOutNA <- rbind(sample1)

    # run the module on fMAF, rSNV, silent, writeLog and na.rm as parameters
    moduleOutNA <- importSNV.TCGA(c("testInNA"), "moduleNA", silent=FALSE, writeLog=FALSE, na.rm=TRUE)
    moduleNA <- read.csv("moduleOutNA", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutNA, moduleNA)
})


test_that("a corrupt input does not lead to corrupted output", {
    # if fMAF file open does not follow the 2.4 or agreed upon version either STOP or DROP file
    expect_error(importSNV.TCGA(c("testIn2.4.1")))
    # if fMAF file open does not have the required columns at the
    #           right column number STOP or DROP file
    expect_error(importSNV.TCGA(c("testInWrongColumns")))
    # if rSNV is provided, check if number of columns in rSNV is worng
    expect_error(importSNV.TCGA(c("testIn")), "testOutINVALIDcolumns")
    # if rSNV is provided, check when header labels are wrong
    expect_error(importSNV.TCGA(c("testIn")), "testOutINVALIDheaders")
})

test_that("silent and write works as intended", {
    # if silent=TRUE, check if output to console is supressed
    testF <- tempfile()
    capture.output(importSNV.TCGA(c("testIn"), silent = TRUE), file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))

    # if writeLog=TRUE, check if log file is created

    # if writeLog=TRUE, check if logs are written to it correctly

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
