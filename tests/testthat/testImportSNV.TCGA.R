# testImportSNV.TCGA.R

# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak
# ==== END SETUP AND PREPARE ===================================================

#context("Tools for importing MAF files into a rSNV file")
options(stringsAsFactors=FALSE)

test_that("parameter errors are correctly handled", {
    # fMAF is not provided
    expect_error(importSNV.TCGA(), "argument \"fMAF\" is missing, with no default")
    # fMAF is NULL/empty vector
    expect_error(importSNV.TCGA(fMAF = NULL), "parameter fMAF cannot take value \"NULL\": should be a vector/string of MAF file name(s)")
    # fMAF input vector contains non-string values
    expect_error(importSNV.TCGA(fMAF = 123), "parameter fMAF cannot be numeric: should be a vector/string of MAF file name(s)")

    # if rSNV parameter is provided, but is NULL
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = NULL), "parameter rSNV cannot take value \"NULL\": should be absent or existing rSNV file name")
    # if rSNV parameter is provided, but is non-string value
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = 123), "parameter rSNV cannot be numeric: should be absent or existing rSNV file name")
    # if rSNV is given vector input
    expect_error(importSNV.TCGA(fMAF = testMAF, rSNV = c("abc", "xyz")), "parameter rSNV cannot take more than one file as input")

    # if silent is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,silent = 1), "silent cannot be numeric : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,silent = NULL), "silent cannot be NULL : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,silent = c(TRUE, TRUE)), "silent cannot have more than 1 values : should be logical TRUE/FALSE")

    # writeLog is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,writeLog = 1), "write.log cannot be numeric : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,writeLog = NULL), "write.log cannot be NULL : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,writeLog = c(TRUE, TRUE)), "write.log cannot have more than 1 values : should be logical TRUE/FALSE")

    # if na.rm is not given logical input
    expect_error(importSNV.TCGA(fMAF = testMAF,na.rm = 1), "na.rm cannot be numeric : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,na.rm = NULL), "na.rm cannot be NULL : should be logical TRUE/FALSE")
    expert_error(importSNV.TCGA(fMAF = testMAF,na.rm = c(TRUE, TRUE)), "na.rm cannot have more than 1 values : should be logical TRUE/FALSE")

})

test_that("a sane input gives an expected output with differnt input parameters", {

    # MAKE SMALL INPUT (testIn) ----------------------------------------------------
    testIn <- tempfile()
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input1 <- c("TMEM201" , 199953 , "WUGSC" , "GRCh38" , "chr1" , 9607661 , 9607661 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a", "c.1265C>T" , "p.Pro422Leu", "p.P422L" , "ENST00000340381")
    input2 <- c("PTCHD2" , 57540 , "WUGSC" , "GRCh38" , "chr1" , 11516130 , 11516130 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.1718C>T" , "p.Ala573Val" , "p.A573V" , "ENST00000294484")
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueIn <- rbind(input0,input1,input2,input3)
    write.table(trueIn,  file=testIn, sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # MAKE SMALL OUTPUT (trueOut) --------------------------------------------------
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")
    sample1 <-c("TMEM201", "chr1", 9607661, 9607661, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("PTCHD2", "chr1", 11516130, 11516130, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA.QQ.A8VF.01A.11D.A37C.09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")

    trueOut <- rbind(header,sample1,sample2,sample3)

    # run the module on only fMAF as parameter (moduleOut)
    moduleOut <- importSNV.TCGA(testIn)
    # check if manual output = module output
    module <- read.table(moduleOut, sep = "\t", stringsAsFactors = FALSE)

    # TEST fMAF parameter ==========================================================
    all.equal(trueOut, module)
    #===============================================================================


    # test input ----------------------------------------------------------------
    #  have an input file, rows too big to import manually here (testInPre)
    testInPre <- tempfile()
    module2 <- tempfile()
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueInPre <- rbind(input0,input3)

    write.table(trueInPre,  file=testInPre, sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # test rSNV input file
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")
    snvInput1 <-c("TMEM201", "chr1", 9607661, 9607661, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    snvInput2 <- c("PTCHD2", "chr1", 11516130, 11516130, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA.QQ.A8VF.01A.11D.A37C.09", "2109690d-e8b1-4b5c-9a39-21af273324cf")

    module <- rbind(header,snvInput1, snvInput2)

    write.table(module,  file=module2, sep="\t", append=FALSE)

    # test output for rSNV as input
    trueOutPre <- trueOut

    # TEST fMAF, rSNV parameters ====================================================
    # run the module on fMAF and rSNV as parameters
    moduleOutPre1 <- importSNV.TCGA(testInPre, module02)
    modulePre1 <- read.table(moduleOutPre1, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre1)
    #===============================================================================

    # TEST fMAF, rSNV, silent parameters ============================================
    # run the module on fMAF and rSNV and silent as parameters
    moduleOutPre2 <- importSNV.TCGA(testInPre, module02, silent=TRUE)
    modulePre2 <- read.table(moduleOutPre2, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre2)
    #===============================================================================

    # TEST fMAF, rSNV, silent, wrtieLog parameters =================================
    # run the module on fMAF, rSNV, silent and writeLog as parameters
    moduleOutPre3 <- importSNV.TCGA(testInPre, module02, silent=FALSE, writeLog=FALSE)
    modulePre3 <- read.table(moduleOutPre3, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre3)
    #===============================================================================

    expect_true(file.remove(testInPre))
    expect_true(file.remove(module2))

    # test input with rSNV as input and NA values -------------------------------------------
    # have an input file, rows too big to import manually here (testInNA)
    testInNA <- tempfile()
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueInNA <- rbind(input0,input3)

    write.table(trueInNA,  file=testInNA, sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # test output for rSNV as input and NA values
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")

    trueOutNA <- rbind(header)

    # TEST fMAF, rSNV, silent, wrtieLog parameters, na.rm ==========================
    # run the module on fMAF, rSNV, silent, writeLog and na.rm as parameters
    moduleOutNA <- importSNV.TCGA(testInNA, silent=FALSE, writeLog=FALSE, na.rm=TRUE)
    moduleNA <- read.table(moduleOutNA, sep = "\t", stringsAsFactors = FALSE)
    # add to prexisting rSNV file
    all.equal(trueOutNA, moduleNA)
    #===============================================================================
    expect_true(file.remove(testInNA))
})


test_that("a corrupt input does not lead to corrupted output", {
    # if fMAF file open does not follow the 2.4 or agreed upon version either STOP or DROP file
    # expect_error(importSNV.TCGA(c("testIn2.4.1")))

    # if fMAF file open does not have the required columns at the
    #           right column number STOP or DROP file
})

test_that("silent and writeLog works as intended", {
    # if silent=TRUE, check if output to console is supressed
    testF <- tempfile()
    capture.output(importSNV.TCGA(testIn, silent = TRUE), file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))
    expect_true(file.remove(testIn))
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
