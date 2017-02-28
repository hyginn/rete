# testFastMap.R
#
#
context("fastMap functions")

test_that("fastMap maps an ID properly", {
    # Generate a small hash table data set
    fastMapUniProt <<- new.env(hash = TRUE)
    fastMapENSP <<- new.env(hash = TRUE)
    fastMapUniProt[["P04217"]] <- "A1BG"
    fastMapUniProt[["Q9NQ94"]] <- "A1CF"
    fastMapUniProt[["P01023"]] <- "A2M"
    fastMapENSP[["ENSG00000268895"]] <- "A1BG-AS1"
    fastMapENSP[["ENSG00000245105"]] <- "A2M-AS1"

    # Single lookup
    expect_equal(fastMap("P04217"), "A1BG")
    expect_equal(fastMap("ENSG00000268895", type = "ENSP"), "A1BG-AS1")

    # Vectorized lookup
    expect_identical(fastMap(c("Q9NQ94", "P01023")), c("A1CF", "A2M"))

    # ID not found
    expect_warning(x <- fastMap(c("P1234", "Q321")))
    expect_identical(x, c("P1234", "Q321"))
})

test_that("fastMapUpdate modifies a fastMap hash table correctly", {
    fastMapUniProt <<- new.env(hash = TRUE)

    # Insert
    fastMapUpdate(fastMapUniProt, "P04217", "A1BG")
    expect_equal(fastMapUniProt[["P04217"]], "A1BG")

    # Modify
    fastMapUpdate(fastMapUniProt, "P04217", "GB1A")
    expect_equal(fastMapUniProt[["P04217"]], "GB1A")

    # Delete
    fastMapUpdate(fastMapUniProt, "P04217", NULL)
    expect_equal(fastMapUniProt[["P04217"]], NULL)

    # Insert multicase
    fastMapUpdate(fastMapUniProt, "Q9nQ94", "a1cF")
    expect_equal(fastMapUniProt[["Q9nQ94"]], "A1CF")
    expect_equal(fastMapUniProt[["Q9NQ94"]], NULL)

    # Insert with acceptable special symbol (- and _)
    fastMapUpdate(fastMapUniProt, "A3-1B", "Z_A92")
    expect_equal(fastMapUniProt[["A3-1B"]], "Z_A92")

    # Insert invalid key
    expect_error(fastMapUpdate(fastMapUniProt, "P123$", "A1BG"))
    expect_error(fastMapUpdate(fastMapUniProt, "", "A1BG"))
    expect_error(fastMapUpdate(fastMapUniProt, NULL, "A1BG"))

    # Insert invalid value
    expect_error(fastMapUpdate(fastMapUniProt, "P04217", "A1 BG"))
})

test_that("fastMapGenerate imports from an HGNC dataset", {
    uniprot <- new.env(hash = TRUE)
    ensp <- new.env(hash = TRUE)

    uniprot[["P04217"]] <- "A1BG"
    uniprot[["Q9NQ94"]] <- "A1CF"
    uniprot[["P01023"]] <- "A2M"
    ensp[["ENSG00000121410"]] <- "A1BG"
    ensp[["ENSG00000268895"]] <- "A1BG-AS1"
    ensp[["ENSG00000148584"]] <- "A1CF"
    ensp[["ENSG00000175899"]] <- "A2M"
    ensp[["ENSG00000245105"]] <- "A2M-AS1"

    fastMapGenerate("hgnc_subset.txt", type = "UniProt", saveHashTable = FALSE)
    fastMapGenerate("hgnc_subset.txt", type = "ENSP", saveHashTable = FALSE)

    expect_equal(uniprot[["P04217"]], fastMapUniProt[["P04217"]])
    expect_equal(uniprot[["Q9NQ94"]], fastMapUniProt[["Q9NQ94"]])
    expect_equal(uniprot[["P01023"]], fastMapUniProt[["P01023"]])
    expect_equal(ensp[["ENSG00000121410"]], fastMapENSP[["ENSG00000121410"]])
    expect_equal(ensp[["ENSG00000268895"]], fastMapENSP[["ENSG00000268895"]])
    expect_equal(ensp[["ENSG00000175899"]], fastMapENSP[["ENSG00000175899"]])

    # Value not inserted
    expect_equal(ensp[["ENSG0"]], fastMapENSP[["ENSG0"]])

    # Conflicting keys
    expect_warning(fastMapGenerate("hgnc_with_conflict.txt", type = "UniProt",
                                   saveHashTable = FALSE))
    expect_equal(fastMapUniProt[["Q9NQ94"]], "A1BG")
})

test_that("fastMap functions integrate together", {
    # Generate the hash table
    fastMapGenerate("hgnc_subset.txt", type = "UniProt", saveHashTable = FALSE)
    # fastMap and fastMapGenerate interaction
    expected <- c("A1BG", "A1CF", "A2M")
    expect_identical(fastMap(c("P04217", "Q9NQ94", "P01023")), expected)
    # fastMap and fastMapUpdate interaction
    fastMapUpdate(fastMapUniProt, "P04217", NULL)
    fastMapUpdate(fastMapUniProt, "Q9NQ94", "PTEST")
    expect_identical(fastMap("P04217"), "P04217")
    expect_identical(fastMap("Q9NQ94"), "PTEST")
})
