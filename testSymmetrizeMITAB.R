# ==== BEGIN SETUP AND PREPARE =================================================
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
NL <- .PlatformLineBreak()

interactor_A = c("fish1", "cat1", "dog1")
interactor_B = c("fish2", "cat2", "dog2")
method = c("MI:0058(genome based prediction)", "MI:1232(aggregation assay)", "MI:0686(unspecified method)")
user_method = c("MI:0058(genome based prediction)", "MI:1232(aggregation assay)", "MI:0686(unspecified method)", "MI:0949(gdp/gtp exchange assay)")
test_interactions = cbind(interactor_A, interactor_B, method)
test_symmetric_methods = cbind(method, c(FALSE, TRUE, NA, TRUE))
interactor_A_corrupt = c("fish1", NA, "dog1")
interactor_B_corrupt = c(NA, NA, NA)
method_corrupt = c(NA, "MI:1232(aggregation assay)", "MI:0686(unspecified method)", "MI:0949(gdp/gtp exchange assay)")
test_interactions_corruptA = cbind(interactor_A_corrupt, interactor_B, methods)
test_interactions_corruptMethods = cbind(interactor_A, interactor_B, method_corrupt)
test_interactions_corruptMix = cbind(interactor_A_corrupt, interactor_B_corrupt, method_corrupt)

# ==== END SETUP AND PREPARE ===================================================

test_that("parameter errors are correctly handled", {
    # Try each parameter missing in turn
    # Try parameters out of expected bounds
    # Try NULL and length-zero parameters: be sure none lead to an erroneous
    #    one-time execution of any loop

    ## test to see if providing correct parameters will give the correct result
    expect_equal(symmetrize(test_intractions, user_method), nrow(test_interactions) == 4 &
                                                        test_interactions$interactor_A[4] == "cat2" &
                                                        test_interactions$interactor_B[4] == "cat1" &
                                                        test_interactions$method == "MI:1232(aggregation assay)")
    ## test to see if both missing parameters gives an error
    expect_error(symmetrize())
    ## test to see if wrong boolean value in first parameter gives an error
    expect_error(symmetrize(TRUE, methods))
    ## test to see if wrong boolean value for the second parameter gives an error
    expect_error(symmetrize(test_interactions, FALSE))
    ## test to see if wrong boolean values in both parameters gives an error
    expect_error(symmetrize(TRUE, TRUE))
    ## test to see if wrong string value for both parameters gives an error
    expect_error(symmetrize("test_interactions", "methods"))
    ## test to see if wrong string value in first parameter gives an error
    expect_error(symmetrize("test_interactions", methods))
    ## test to see if worng string value in second parameter gives an error
    expect_error(symmetrize(test_interactions, "methods"))
    ## test to see if NULL value in second parameter gives an error
    expect_error(symmetrize(test_interactions, NULL))
    ## test to see if NULL value in first parameter gives an error
    expect_error(symmetrize(NULL, methods))
    ## test to see if NULL values in both parameters gives an error
    expect_error(symmetrize(NULL, NULL))
    ## test to see if missing value for the second parameter gives an error
    expect_error(symmetrize(test_interactions))
    ## test to see if missing value for the first parameter gives an error
    expect_error(symmetrize(methods))
    ## test to see if NULL being the only parameter gives an error
    expect_error(symmetrize(NULL))
    ##  test to see if adding an extra parameter gives an error
    expect_error(symmetrize(test_interactions, methods, extra))
})


test_that("a sane input gives an expected output", {
    # Provide small input data, provide small output data.
    #   (Note that the output needs to be independently verifiably correct.
    #    Don't run an erroneous function, and then test against the erroneous
    #    output you received when you ran your function for the first time.)
    # Try this with important variations of parameters (reasonable, not
    #    combinatorially exhaustive).
    # Cover your code.

    ## test to see if function will correctly add a symmetric edge to the input
    expect_equal(symmetrize(test_intractions, user_method), nrow(test_interactions) == 4 &
                     test_interactions$interactor_A[4] == "cat2" &
                     test_interactions$interactor_B[4] == "cat1" &
                     test_interactions$method == "MI:1232(aggregation assay)")
    ## test to see if function will correctly not add a symmetric edge to the input, when the symmetric edge is not required
    expect_equal(symmetrize(test_intractions[1, ], user_method), nrow(test_interactions) == 1 &
                     test_interactions$interactor_A == "fish1" $
                     test_interactions$interactor_B == "fish2" $
                     test_interactions$method == "MI:0058(genome based prediction)")
    ## test to see if function will correctly not add a symmeric edge if the type of method is not specific enough to determine if
    ## a symmetric edge is required
    expect_equal(symmetrize(test_interactions[3, ], user_method), nrow(test_interactions) == 1 &
                     test_interactions$interactor_A == "dog1" &
                     test_interactions$inderactor_B == "dog2" &
                     test_interactions$method == "MI:0686(unspecified method)")
    ## test to see if the function correctly removes interactions with methods that do not match the requested user methods
    ## also checks to see if the function corretly does not add a symmetric edge in this case
    expect_equal(symmetrize(test_interactions, user_method[1]), nrow(test_interactions) == 1 &
                     test_interactions$interactor_A == "fish1" &
                     test_interactions$inderactor_B == "fish2" &
                     test_interactions$method == "MI:0058(genome based prediction)")
    ## test to see if the function correctly removes interactions with methods that do not match the requested user methods
    ## also checks to see if the function corretly adds a symmetric edge, since it is required in this case
    expect_equal(symmetrize(test_interactions, user_method[2]), nrow(test_interactions) == 2 &
                     test_interactions$interactor_A[1] == "cat1" &
                     test_interactions$inderactor_B[1] == "cat2" &
                     test_interactions$method == "MI:0058(genome based prediction)" &
                     test_interactions$interactor_A[2] == "cat2" &
                     test_interactions$inderactor_B[2] == "cat1" &
                     test_interactions$method == "MI:0058(genome based prediction)")
    ## test to see if the function correctly removes all interactions with methods that do not match the user requested methods
    ## in this case, no interaction methods should match the user requested methods, so all interactions should be removed
    expect_equal(symmetrize(test_interactions, user_method[4]), nrow(test_interactions) == 0)
})


test_that("a corrupt input does not lead to corrupted output", {
    # Try: - spurious characters in numeric columns ("N/A", ...).
    #      - extra tabs at line-end
    #      - comments before headers
    #      - truncated file (incomplete last line)
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop

    ## test to see if an input with a corrupt column in the first parameter, correctly results in an error
    expect_error(symmetrize(test_interactions_corruptA, user_method))
    ## test to see if an input with a corrupt second parameter correctly results in an error
    expect_error(symmetrize(test_interactions, method_corrupt))
    ## test to see if an input with a corrupt methods column in the first parameter, correctly results in an error
    expect_error(symmetrize(test_interactions_corruptMethods, user_method))
    ## test to see if an input with corrupt columns in all of the first parameter correctly results in an error
    expect_error(symmetrize(test_interactions_corruptMix, user_method))
    ## test to see if an input with a corrupt first parameter and second parameter correctly results in an error
    expect_error(symmetrize(test_interactions_corruptMix, method_corrupt))
})



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}
options("rete.logfile" = OLOG)
rm(test_interactions, test_interactions_corruptMix, test_interactions_corruptMethods, test_interactions_corruptA,
   test_symmetric_methods, method, method_corrupt, interactor_B_corrupt, interactor_A_corrupt, interactor_B, interactor_A, user_method)
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
