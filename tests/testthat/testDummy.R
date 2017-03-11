# testDummy.R
#
# Dummy tests as a template for unit testing in test-driven development
# for the rete project.
context("<Dummy tests>")

# for demonstration only - remove!
dummyF <- function(a) {
    if (is.null(a) || length(a) == 0 || !is.numeric(a)) { stop() }
    attr(a, "date") <- Sys.Date()
    return(a * a)
}


test_that("parameter errors are correctly handled", {
    # Try each parameter missing in turn
    # Try parameters out of expected bounds
    # Try NULL and length-zero parameters: be sure none lead to an erroneous
    #    one-time execution of any loop
    expect_equal(as.numeric(dummyF(2)), 4)
    expect_equal(as.numeric(dummyF(-1/3)), 1/9)
    expect_equal(as.numeric(dummyF(-1/0)), Inf)
    expect_error(dummyF())
    expect_error(dummyF(TRUE))
    expect_error(dummyF("2"))
    expect_error(dummyF(NULL))
    expect_error(dummyF(numeric()))
})


test_that("a sane input gives an expected output", {
    # Provide small input data, provide small output data.
    #   (Note that the output needs to be independently verifiably correct.
    #    Don't run an erroneous function, and then test against the erroneous
    #    output you received when you ran your function for the first time.)
    # Try this with important variations of parameters (reasonable, not
    #    combinatorially exhaustive).
    # Cover your code.
    expect_equal(as.numeric(dummyF(2)), 4)
    expect_equal(attr(dummyF(2), "date"), Sys.Date())
    expect_equal(dummyF(dummyF(2)), dummyF(4))
})


test_that("a corrupt input does not lead to corrupted output", {
    # Try: - spurious characters in numeric columns ("N/A", ...).
    #      - extra tabs at line-end
    #      - comments before headers
    #      - truncated file (incomplete last line)
    # Again: make absolutely sure be you never have an erroneous
    #        one-time execution of a loop
    expect_error(dummyF(NA))
})


# [END]
