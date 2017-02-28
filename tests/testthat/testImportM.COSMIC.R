# testImportM.COSMIC.R
#
#
context("import COSMIC mutations")

test_that("importM.COSMIC correctly filters TSV file with NCV data", {

    importM.COSMIC(fName = "testCosmicNCV.tsv",
                   type = "NCV",
                   range <- c(1:2),
                   verbose = FALSE)

    expected <- data.frame("Sample name" = c("012-02-1TD",
                                             "PD4771a",
                                             "GC8_T",
                                             "GC8_T",
                                             "GC8_T",
                                             "CHC121T",
                                             "PD3905a",
                                             "PD4103a",
                                             "PD4109a",
                                             "LFS_MB1"),
                           "ID_SAMPLE" = c(1529895,
                                           1432693,
                                           1645684,
                                           1645684,
                                           1645684,
                                           1652968,
                                           1230729,
                                           1317020,
                                           1317053,
                                           1652860),
                           stringsAsFactors = FALSE)

    expect_equal(read.table("./NCV.txt",
                            sep="\t",
                            header=TRUE,
                            stringsAsFactors = FALSE),
                 expected)
})

test_that("importM.COSMIC correctly filters TSV file with CNV data", {

    importM.COSMIC(fName = "testCosmicCompleteCNA.tsv",
                   type = "CNV",
                   range <- c(1:2),
                   verbose = FALSE)

    expected <- data.frame("CNV_ID" = c(8131112,
                                        7894610,
                                        7969471,
                                        8206785,
                                        8132439,
                                        8206785,
                                        8170609,
                                        8170609,
                                        8071525,
                                        7984510),
                           "ID_GENE" = c(69775,
                                         75764,
                                         66768,
                                         84805,
                                         68055,
                                         84836,
                                         103307,
                                         69786,
                                         69782,
                                         106410),
                           stringsAsFactors = FALSE)

    expect_equal(read.table("./CNV.txt",
                            sep="\t",
                            header=TRUE,
                            stringsAsFactors = FALSE),
                 expected)
})

test_that("importM.COSMIC returns error when file not found", {

    expect_error(
        importM.COSMIC(fName = "testCosmicNCV2.tsv",
                       type = "NCV",
                       range <- c(1:2),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when type not 'NCV' or 'CNV'", {

    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV2",
                       range <- c(1:2),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when zero index in range", {

    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- c(1:2, 0),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when less than zero index in range", {

    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- c(1:2, -1),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when index > ncols in fName table", {

    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- c(19),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when fName is not a valid string", {
    expect_error(
        importM.COSMIC(fName = 0,
                       type = "NCV",
                       range <- c(1:2),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when type is not a valid string", {
    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = 0,
                       range <- c(1:2),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when range is not a valid vector", {
    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- 0,
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when index is not a valid numeric", {
    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- c(1:2, "3"),
                       verbose = FALSE))
})

test_that("importM.COSMIC returns error when verbose is not a valid logical", {
    expect_error(
        importM.COSMIC(fName = "testCosmicNCV.tsv",
                       type = "NCV",
                       range <- c(1:2),
                       verbose = "FALSE"))
})

test_that("importM.COSMIC outputs entire table when range not given", {

    importM.COSMIC(fName = "testSmallCosmicNCV.tsv",
                   type = "NCV",
                   verbose = FALSE)

    expected <- data.frame("Sample name" =
                               c("012-02-1TD",
                                 "PD4771a"
                               ),
                           "ID_SAMPLE" =
                               c(1529895,
                                 1432693
                               ),
                           "ID_NCV" =
                               c("COSN159757",
                                 "COSN161488"
                               ),
                           "zygosity" =
                               c("Heterozygous",
                                 "Heterozygous"
                               ),
                           "GRCh" =
                               c(38,
                                 38
                               ),
                           "genome position" =
                               c("5:155388882-155388882",
                                 "20:32889328-32889328"
                               ),
                           "Mutation somatic status" =
                               c("Confirmed somatic variant",
                                 "Confirmed somatic variant"
                               ),
                           "WT_SEQ" =
                               c("C",
                                 "C"
                               ),
                           "MUT_SEQ" =
                               c(TRUE,
                                 TRUE
                               ),
                           "SNP" =
                               c("n",
                                 "n"
                               ),
                           "FATHMM_MKL_NON_CODING_SCORE" =
                               c(0.09499,
                                 0.97855
                               ),
                           "FATHMM_MKL_NON_CODING_GROUPS" =
                               c("A",
                                 "AD"
                               ),
                           "FATHMM_MKL_CODING_SCORE" =
                               c(0.00874,
                                 0.79666
                               ),
                           "FATHMM_MKL_CODING_GROUPS" =
                               c("AEFI",
                                 "AEFDGI"
                               ),
                           "Whole_Genome_Reseq" =
                               c("y",
                                 "n"
                               ),
                           "Whole_Exome" =
                               c("n",
                                 "y"
                               ),
                           "ID_STUDY" =
                               c(340,
                                 351
                               ),
                           "PUBMED_PMID" =
                               c(21642962,
                                 21995386
                               ),
                           stringsAsFactors = FALSE)

    expect_equal(read.table("./NCV.txt",
                            sep="\t",
                            header=TRUE,
                            stringsAsFactors = FALSE),
                 expected)
})

# [END]
