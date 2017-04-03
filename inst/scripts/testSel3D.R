# testSel3D.R
#
#
# Structure your tests according to the three principles expressed below:
#    test_that("parameter errors are correctly handled", { ... })
#    test_that("a sane input gives an expected output", { ... })
#    test_that("a corrupt input does not lead to corrupted output", { ...})


# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()
# ==== END SETUP AND PREPARE ===================================================


someResults <- "MLLLARCLLLVLVSSLLVCSGLACGPGRGFGKRRHPKKLTPLAYKQFIPNVAEKTLGASGRYEGKISRNSERFKELTPNYNPDIIFKDEENTGADRLMTQRCKDKLNALAISVMNQWPGVKLRVTEGWDEDGHHSEESLHYEGRAVDITTSDRDRSKYGMLARLAVEAGFDWVYYESKAHIHCSVKAENSVAAKSGGCFPGSATVHLEQGGTKLVKDLSPGDRVLAADDQGRLLYSDFLTFLDRDDGAKKVFYVIETREPRERLLLTAAHLLFVAPHNDSATGEPEASSGSGPPSGGALGPRALFASRVRPGQRVYVVAERDGDRRLLPAAVHSVTLSEEAAGAYAPLTAQGTILINRVLASCYAVIEEHSWAHRAFAPFRLAHALLAALAPARTDRGGDSGGGDRGGGGGRVALTAPGAADAPGAGATAGIHWYSQLLYQIGTWLLDSEALHPLGMAVKSS"
manyResults <- "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
test_that("parameter errors are correctly handled", {
  # Try missing input sequence
  # Try extraneous parameters
  expect_error(SEL3D())
  expect_error(SEL3D(""))
  expect_error(SEL3D("AGT", 3))
  expect_error(SEL3D("AGT", silent=1), "Mode error: argument requires logical input")
})


test_that("a sane input gives an expected output", {
  # Note: All of the expected results need to be verified by BLAST to PDB
  # Try single codon
  # Try an input which has no results (verify with BLAST to PDB)
  # Try small input with a small result
  # Try a large input with many resulting matches (verify with BLAST to PDB)

  expect_error(SEL3D("AGT"))
  expect_error(SEL3D("AGTCTTATTGTTTTT"))
  #although the above returns a result on BLAST, the module specification dictates
  #that we shouldn't return results for amino acid sequences who's length is <50

  #Manually check if top matches on BLAST are present in the output for the right positions
  result <- SEL3D(someResults, convert=FALSE)
  expect_equal(length(result$pdb.id), 1)
  expect_true("3M1N_A" %in% result$pdb.id)
  #check q.start
  expect_equal(result[result$pdb.id == "3M1N_A", 8], 25)
  #check q.end
  expect_equal(result[result$pdb.id == "3M1N_A", 9], 197)
  #check s.start
  expect_equal(result[result$pdb.id == "3M1N_A", 10], 3)
  #check s.end
  expect_equal(result[result$pdb.id == "3M1N_A", 11], 175)

  #check that the silent param works
  expect_message(SEL3D(someResults, convert=FALSE, silent=FALSE))


  #Manually check if top matches on BLAST are present in the output for the right positions
  result <- SEL3D(manyResults, convert=FALSE)
  expect_equal(length(result$pdb.id), 4)
  expect_true("4QO1_B" %in% result$pdb.id)
  expect_equal(result[result$pdb.id == "4QO1_B", 8], 92)
  expect_equal(result[result$pdb.id == "4QO1_B", 9], 312)
  expect_equal(result[result$pdb.id == "4QO1_B", 10], 1)
  expect_equal(result[result$pdb.id == "4QO1_B", 11], 221)
  expect_equal(result[result$pdb.id == "4QO1_B", 17], "Homo sapiens")

  expect_true("5HPD_A" %in% result$pdb.id)
  expect_equal(result[result$pdb.id == "5HPD_A", 8], 2)
  expect_equal(result[result$pdb.id == "5HPD_A", 9], 61)
  expect_equal(result[result$pdb.id == "5HPD_A", 10], 99)
  expect_equal(result[result$pdb.id == "5HPD_A", 11], 158)
  expect_equal(result[result$pdb.id == "5HPD_A", 17], "Homo sapiens")

  expect_true("1OLG_A" %in% result$pdb.id)
  expect_equal(result[result$pdb.id == "1OLG_A", 8], 319)
  expect_equal(result[result$pdb.id == "1OLG_A", 9], 360)
  expect_equal(result[result$pdb.id == "1OLG_A", 10], 1)
  expect_equal(result[result$pdb.id == "1OLG_A", 11], 42)
  expect_equal(result[result$pdb.id == "1OLG_A", 17], "Homo sapiens")

  expect_true("1DT7_X" %in% result$pdb.id)
  expect_equal(result[result$pdb.id == "1DT7_X", 8], 367)
  expect_equal(result[result$pdb.id == "1DT7_X", 9], 388)
  expect_equal(result[result$pdb.id == "1DT7_X", 10], 1)
  expect_equal(result[result$pdb.id == "1DT7_X", 11], 22)
  expect_equal(result[result$pdb.id == "1DT7_X", 17], "Homo sapiens")

})


test_that("a corrupt input does not lead to corrupted output", {
  # Try illegal characters in input sequence
  # Try illegal length of input sequence
  # Try empty string as input sequence
  # Try unexpected type of input sequence
  # These should all exit immediately without executing any of the actual algo.
  expect_error(SEL3D("AAAA"))
  expect_error(SEL3D("AGTCXT"))
  expect_error(SEL3D(""))
  expect_error(SEL3D(c("A", "G", "T")))
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
