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

    # make small input (testIn) ----------------------------------------------------
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input1 <- c("TMEM201" , 199953 , "WUGSC" , "GRCh38" , "chr1" , 9607661 , 9607661 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a", "c.1265C>T" , "p.Pro422Leu", "p.P422L" , "ENST00000340381")
    input2 <- c("PTCHD2" , 57540 , "WUGSC" , "GRCh38" , "chr1" , 11516130 , 11516130 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09" , "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.1718C>T" , "p.Ala573Val" , "p.A573V" , "ENST00000294484")
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueIn <- rbind(input0,input1,input2,input3)
    write.table(trueIn,  file="testIn", sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # make small output for that input manually (testOut)
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")
    sample1 <-c("TMEM201", "chr1", 9607661, 9607661, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("PTCHD2", "chr1", 11516130, 11516130, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA.QQ.A8VF.01A.11D.A37C.09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")

    trueOut <- rbind(header,sample1,sample2,sample3)

    # run the module on only fMAF as parameter (moduleOut)
    moduleOut <- importSNV.TCGA(c("testIn"))
    # check if manual output = module output
    module <- read.csv(moduleOut, sep = "\t")
    all.equal(trueOut, module)

    # test input with rSNV as input ---------------------------------------------
    #  have an input file, rows too big to import manually here (testInPre)
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueInPre <- rbind(input0,input3)

    write.table(trueInPre,  file="testInPre", sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # test output for rSNV as input
    trueOutPre <- trueOut

    # run the module on fMAF and rSNV as parameters
    moduleOutPre1 <- importSNV.TCGA(c("testInPre"), "module02")
    modulePre1 <- read.csv(moduleOutPre1, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre1)

    # run the module on fMAF and rSNV and silent as parameters
    moduleOutPre2 <- importSNV.TCGA(c("testInPre"), "module02", silent=TRUE)
    modulePre2 <- read.csv(moduleOutPre2, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre2)

    # run the module on fMAF, rSNV, silent and writeLog as parameters
    moduleOutPre3 <- importSNV.TCGA(c("testInPre"), "module02", silent=FALSE, writeLog=TRUE)
    modulePre3 <- read.csv(moduleOutPre3, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre3)

    # test input with rSNV as input and NA values -------------------------------------------
    #  have an input file, rows too big to import manually here (testInNA)
    input0 <- c("Hugo_Symbol" , "Entrez_Gene_Id" , "Center" , "NCBI_Build" , "Chromosome" , "Start_Position" , "End_Position" , "Strand" , "Variant_Classification" , "Variant_Type" , "Reference_Allele" , "Tumor_Seq_Allele1" , "Tumor_Seq_Allele2" , "dbSNP_RS" , "dbSNP_Val_Status" , "Tumor_Sample_Barcode" , "Matched_Norm_Sample_Barcode" , "Match_Norm_Seq_Allele1" , "Match_Norm_Seq_Allele2" , "Tumor_Validation_Allele1" , "Tumor_Validation_Allele2" , "Match_Norm_Validation_Allele1" , "Match_Norm_Validation_Allele2" , "Verification_Status" , "Validation_Status" , "Mutation_Status" )
    input3 <- c("VWA5B1" , 127731 , "WUGSC" , "GRCh38" , "chr1" , 20319393 , 20319393 , "+" , "Missense_Mutation" , "SNP" , "C" , "C" , "T" , "novel" , "TCGA-QQ-A8VF-01A-11D-A37C-09", "TCGA-QQ-A8VF-10A-01D-A37F-09" , "Somatic" , "Illumina" , "HiSeq" , "2000" , "2109690d-e8b1-4b5c-9a39-21af273324cf" , "fe2daf93-56f1-4236-9b58-bf18f0e47d3a" , "c.853C>T" , "p.Pro285Ser" , "p.P285S" , "ENST00000375079")

    trueInNA <- rbind(input0,input3)

    write.table(trueInNA,  file="testInNA", sep="\t", append=FALSE, row.names =
                    FALSE, col.names = FALSE, quote = FALSE)

    # test output for rSNV as input and NA values
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")

    trueOutNA <- rbind(header)

    # run the module on fMAF, rSNV, silent, writeLog and na.rm as parameters
    moduleOutNA <- importSNV.TCGA(c("testInNA"), silent=FALSE, writeLog=FALSE, na.rm=TRUE)
    moduleNA <- read.csv(moduleOutNA, sep = "\t")
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

})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================
# [END]
