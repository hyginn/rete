# testImportSNV.TCGA.R

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak
# ==== END SETUP AND PREPARE ===================================================

#context("Tools for importing MAF files into a rSNV file")
options(stringsAsFactors=FALSE)
testIN <- "importSNV_input_test_file.maf"
testOUT <- "importSNV_output_test_file"
testINcorr <- "importSNV_input_corr_file.maf"

test_that("parameter errors are correctly handled", {
    testMAF <- tempfile()
    tOut <- tempfile()
    trueIn <- data.frame()
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input1 <- c("TMEM201" , 199953 , "WUGSC" , "GRCh38" , "chr1" , 9607661 , 9607661 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a", "c.1265C>T" , "p.Pro422Leu", "p.P422L" , "ENST00000340381")
    input2 <- c("PTCHD2" , 57540 , "WUGSC" , "GRCh38" , "chr1" , 11516130 , 11516130 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.1718C>T" , "p.Ala573Val" , "p.A573V" , "ENST00000294484")
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueIn <- rbind(input0, input1, input2, input3)

    write.table(trueIn,  file=testMAF, sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)
    # fMAF is not provided
    expect_error(importSNV.TCGA(), "argument \"fMAF\" is missing, with no default")
    # fMAF is NULL/empty vector
    expect_error(importSNV.TCGA(fMAF = NULL, tOut), ".checkArgs> \"fMAF\" mode error: argument has mode \"NULL\" but function expects mode \"character\".")
    # fMAF input vector contains non-string values
    expect_error(importSNV.TCGA(fMAF = 123, tOut), ".checkArgs> \"fMAF\" mode error: argument has mode \"numeric\" but function expects mode \"character\".")

    # rSNV is missing
    expect_error(importSNV.TCGA(fMAF = testMAF), "argument \"rSNV\" is missing, with no default")
    # if rSNV parameter is provided, but is NULL
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = NULL), ".checkArgs> \"rSNV\" mode error: argument has mode \"NULL\" but function expects mode \"character\".")
    # if rSNV parameter is provided, but is non-string value
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, rSNV = 123), ".checkArgs> \"rSNV\" mode error: argument has mode \"numeric\" but function expects mode \"character\".")
    # if rSNV is given vector input
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, rSNV = c("abc", "xyz")), ".checkArgs> \"rSNV\" length error: argument has length 2 but function expects length 1.")

    # if silent is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, silent = 1), ".checkArgs> \"silent\" mode error: argument has mode \"numeric\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, silent = NULL), ".checkArgs> \"silent\" mode error: argument has mode \"NULL\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, silent = c(TRUE, TRUE)), ".checkArgs> \"silent\" length error: argument has length 2 but function expects length 1.")

    # writeLog is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, writeLog = 1), ".checkArgs> \"writeLog\" mode error: argument has mode \"numeric\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, writeLog = NULL), ".checkArgs> \"writeLog\" mode error: argument has mode \"NULL\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, writeLog = c(TRUE, TRUE)), ".checkArgs> \"writeLog\" length error: argument has length 2 but function expects length 1.")

    # if na.rm is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, na.rm = 1), ".checkArgs> \"na.rm\" mode error: argument has mode \"numeric\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, na.rm = NULL), ".checkArgs> \"na.rm\" mode error: argument has mode \"NULL\" but function expects mode \"logical\".")
    expect_error(importSNV.TCGA(fMAF = testMAF, tOut, na.rm = c(TRUE, TRUE)), ".checkArgs> \"na.rm\" length error: argument has length 2 but function expects length 1.")
    # ToDo - confirm that no log file was written.
    expect_true(file.remove(testMAF))
})

test_that("a sane input gives an expected output with different input parameters", {

    tOut <- tempfile()
    tOut1 <- tempfile()
    tOut2 <- tempfile()
    trueOut <- data.frame()

    # MAKE SMALL OUTPUT (trueOut) --------------------------------------------------
    trueOut <- read.table(testOUT, sep = "\t", stringsAsFactors = FALSE)

    # run the module on only fMAF and rSNV as parameters (moduleOut)
    importSNV.TCGA(testIN, tOut, writeLog = FALSE)
    # check if manual output = module output
    module <- read.table(tOut, sep = "\t", stringsAsFactors = FALSE)

    # TEST fMAF parameter ==========================================================
    expect_equal(trueOut, module)
    #===============================================================================

    # TEST fMAF, rSNV, silent parameters ============================================
    # run the module on fMAF and rSNV and silent as parameters
    importSNV.TCGA(testIN, tOut1, silent=TRUE, writeLog = FALSE)
    modulePre2 <- read.table(tOut1, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    expect_equal(trueOut, modulePre2)
    #===============================================================================

    # TEST fMAF, rSNV, silent, writeLog parameters =================================
    # run the module on fMAF, rSNV, silent and writeLog as parameters
    importSNV.TCGA(testIN, tOut2, silent=FALSE, writeLog=FALSE)
    modulePre3 <- read.table(tOut2, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    expect_equal(trueOut, modulePre3)
    #===============================================================================

    # test input with rSNV as input and NA values -------------------------------------------
    tOut3 <- tempfile()

    # TEST fMAF, rSNV, silent, writeLog parameters, na.rm ==========================
    # run the module on fMAF, rSNV, silent, writeLog and na.rm as parameters
    importSNV.TCGA(testIN, tOut3, silent=FALSE, writeLog=FALSE, na.rm=TRUE)
    moduleNA <- read.table(tOut3, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    expect_equal(trueOut, moduleNA)
})


test_that("a corrupt input does not lead to corrupted output", {
    # checking for version is not necessary : can handle all versions (2.4.1(latest) and below)
    # ToDo: test incomplete line
    tOut <- tempfile()
    expect_warning(importSNV.TCGA(testINcorr, tOut), "incomplete final line found by readTableHeader on \"importSNV_input_corr_file.maf\"")
})

test_that("silent and writeLog work as intended", {
    # if silent=TRUE, check if output to console is supressed
    tOut <- tempfile()
    testF <- tempfile()
    capture.output(importSNV.TCGA(testIN, tOut, silent = TRUE, writeLog = FALSE), file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))

    # check writeLog works as intended
    importSNV.TCGA(testIN, tOut, writeLog = TRUE)
    expect_equal(readLines(unlist(getOption("rete.logfile")))[1], "event | title  | importSNV.TCGA")
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
