# testImportSNV.TCGA.R

context("Tools for importing MAF files into a rSNV file")

test_that("parameter errors are correctly handled", {
    # fMAF is not provided
    expect_error(importSNV())
    # fMAF is NULL/empty vector
    expect_error(importSNV(fMAF = NULL))
    # fMAF input vector contains non-string values
    expect_error(importSNV(fMAF = 123))
    # file in fMAF could not be opened/does not exist
    expect_error(importSNV(fMAF = c("nosuchfileMAF")))

    # if rSNV parameter is provided, rSNV file does not exist
    expect_error(importSNV(fMAF = testMAF, rSNV = "nosuchfileSNV"))
    # if rSNV parameter is provided, but is NULL
    expect_error(importSNV(fMAF = testMAF, rSNV = NULL))
    # if rSNV parameter is provided, but is non-string value
    expect_error(importSNV(fMAF = testMAF, rSNV = 123))
    # if rSNV is given vector input
    expect_error(importSNV(fMAF = testMAF, rSNV = c("abc", "xyz")))
    # if rSNV is not provided, new file could not be created (file name conflict)
    expect_error(importSNV(fMAF = testMAF))

    # if na.rm is not given logical input
    expect_error(importSNV(na.rm = 1))
    expert_error(importSNV(na.rm = NULL))
    expert_error(importSNV(na.rm = logical()))
    expert_error(importSNV(na.rm = c(TRUE, TRUE)))

    # if silent is not given logical input
    expect_error(importSNV(silent = 1))
    expert_error(importSNV(silent = NULL))
    expert_error(importSNV(silent = logical()))
    expert_error(importSNV(silent = c(TRUE, TRUE)))

    # write.log is not given logical input
    expect_error(importSNV(writeLog = 1))
    expert_error(importSNV(writeLog = NULL))
    expert_error(importSNV(wrtieLog = logical()))
    expert_error(importSNV(writeLog = c(TRUE, TRUE)))
})

test_that("a sane input gives an expected output with differnt input parameters", {

    # make small input (testIn)
    testIn <- "~/Desktop/rete/testIn"
    # make small output for that input manually (testOut)
    write("# metadata \n# 1 \n# 2  \n# 3 \n# 4 \n# 5 \n# metadata END\n", "testOut")
    header <- c("sym", "chr", "start", "end", "strand", "class", "type", "aRef", "a1", "a2", "TCGA", "UUID")
    write(header, "testOut", ncolumns = 12, sep = "\t", append = TRUE)
    sample1 <-c("TMEM201", "chr1", 9607661, 9607661, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("PTCHD2", "chr1", 11516130, 11516130, "+", "Missense_Mutation", "SNP", "C", "C", "T", "TCGA.QQ.A8VF.01A.11D.A37C.09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample4 <- c("TXLNA", "chr1",	32192719,	32192719,	"+",	"Silent",	"SNP",	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample5 <- c("NFYC", "chr1",	40770537,	40770537,	"+",	"Silent",	"SNP",	"C",	"C",	"A", "TCGA-QQ-A8VF-01A-11D-A37C-09",	"2109690d-e8b1-4b5c-9a39-21af273324cf")
    write(sample1, "testOut", ncolumns=12, sep="\t", append=TRUE)
    write(sample2, "testOut", ncolumns=12, sep="\t", append=TRUE)
    write(sample3, "testOut", ncolumns=12, sep="\t", append=TRUE)
    write(sample4, "testOut", ncolumns=12, sep="\t", append=TRUE)
    write(sample5, "testOut", ncolumns=12, sep="\t", append=TRUE)
    trueOut <- read.csv("testOut", skip=7, sep = "\t")

    # run the module on only fMAF as parameter (moduleOut)
    moduleOut <- importSNV(c("testIn"))
    # check if manual output = module output
    # skip 7 for metadata lines just for now
    module <- read.csv("moduleOut", skip=7, sep = "\t")
    all.equal(trueOut, module)

    # test input with rSNV as input
    testInPre <- "~/Desktop/rete/testInPre"
    # test output for rSNV as input
    sample1 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample2 <- c("TXLNA", "chr1",	32192719,	32192719,	"+",	"Silent",	"SNP",	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    sample3 <- c("NFYC", "chr1",	40770537,	40770537,	"+",	"Silent",	"SNP",	"C",	"C",	"A", "TCGA-QQ-A8VF-01A-11D-A37C-09",	"2109690d-e8b1-4b5c-9a39-21af273324cf")
    write(sample1, "testOutPre", ncolumns=12, sep="\t", append=TRUE)
    write(sample2, "testOutPre", ncolumns=12, sep="\t", append=TRUE)
    write(sample3, "testOutPre", ncolumns=12, sep="\t", append=TRUE)
    trueOutPre <- read.csv("testOutPre", skip=7, sep = "\t")

    # run the module on fMAF and rSNV as parameters
    moduleOutPre <- importSNV(c("testInPre"), "module02")
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # run the module on fMAF and rSNV and silent as parameters
    moduleOutPre <- importSNV(c("testInPre"), "module02", silent=TRUE)
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # run the module on fMAF, rSNV, silent and writeLog as parameters
    moduleOutPre <- importSNV(c("testInPre"), "module02", silent=TRUE, writeLog=TRUE)
    modulePre <- read.csv("moduleOutPre", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutPre, modulePre)

    # test input with rSNV as input and NA values
    testInNA <- "~/Desktop/rete/testInNA"
    # test output for rSNV as input and NA values
    sample1 <- c("VWA5B1", "chr1",	20319393,	20319393,	"+",	"Missense_Mutation",	"SNP", 	"C",	"C",	"T", "TCGA-QQ-A8VF-01A-11D-A37C-09", "2109690d-e8b1-4b5c-9a39-21af273324cf")
    write(sample1, "testOutNA", ncolumns=12, sep="\t", append=TRUE)
    trueOutNA <- read.csv("testOutNA", skip=7, sep = "\t")

    # run the module on fMAF, rSNV, silent, writeLog and na.rm as parameters
    moduleOutNA <- importSNV(c("testInNA"), "moduleNA", silent=TRUE, writeLog=TRUE, na.rm=TRUE)
    moduleNA <- read.csv("moduleOutNA", skip=7, sep = "\t")
    # add to prexisting rSNV file
    all.equal(trueOutNA, moduleNA)
})


test_that("a corrupt input does not lead to corrupted output", {
    # if fMAF file open does not follow the 2.4 or agreed upon version either STOP or DROP file
    expect_error(importSNV(c("testIn2.4.1")))
    # if fMAF file open does not have the required columns at the
    #           right column number STOP or DROP file
    expect_error(importSNV(c("testInWrongColumns")))
    # if rSNV is provided, check if number of columns in rSNV is worng
    expect_error(importSNV(c("testIn")), "testOutINVALIDcolumns")
    # if rSNV is provided, check when header labels are wrong
    expect_error(importSNV(c("testIn")), "testOutINVALIDheaders")
})

test_that("silent and write works as intended", {
    # if silent=TRUE, check if output to console is supressed
    testF <- tempfile()
    capture.output(Do.A.Thing(silent = TRUE), file = testF)
    expect_equal(length(readLines(testF)), 0)
    expect_true(file.remove(testF))

    # if writeLog=TRUE, check if log file is created

    # if writeLog=TRUE, check if logs are written to it correctly

})


test_that("test that columns contains valid values", {
    # check if variant_class contains valid values from known list
    ### will finalize these tests once import module is implemented ###

    classColumn <- addColumn("testOut")
    classV <- c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation",
                "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation",
                "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR1" , "Intron", "RNA", "Targeted_Region")
    out.C <- TRUE
    i <- 1
    while(out.C == TRUE & i <= length(classColumn)){
        out.C <- is.element(classColumn[i], classV)
        i <- i + 1
    }
    if (i == lenght(classColumn)) {
        print("All values Valid")
    } else {
        print("Invalid values in class column!")
    }

    # check if vairant_type contains valid values from known list
    typeColumn <- addColumn("testOut")
    typeV <- c("SNP", "DNP", "TNP", "ONP", "INS", "DEL", "Consolidated2")
    out.V <- TRUE
    i <- 1
    while(out.V == TRUE & i <= length(classColumn)){
        out.V <- is.element(typeColumn[i], typeV)
        i <- i + 1
    }
    if (i == lenght(typeColumn)) {
        print("All values Valid")
    } else {
        print("Invalid values in type column!")
    }

    # check chromosome value
    chrColumn <- addColumn("testOut")
    typeChr <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13",
                 "chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X","Y","Z")
    out.Chr <- TRUE
    i <- 1
    while(out.Chr == TRUE & i <= length(chrColumn)){
        out.Chr <- is.element(typeColumn[i], typeChr)
        i <- i + 1
    }
    if (i == lenght(chrColumn)) {
        print("All values Valid")
    } else {
        print("Invalid values in chr column!")
    }

    # check start position and end are integers
})

# [END]
