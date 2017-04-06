# testIsHGNCsymbol.R

context("test the closure isHGNCsymbol()")

test_that("expected input is correctly handled", {
    expect_equal(isHGNCsymbol(), logical())
    expect_equal(isHGNCsymbol(NULL), logical())
    expect_false(isHGNCsymbol(""))
    expect_false(isHGNCsymbol(0))
    expect_true(isHGNCsymbol("A1BG"))                    # First in table
    expect_true(isHGNCsymbol("a1bg"))                    # Case insensitive
    expect_true(isHGNCsymbol("ZZZ3"))                    # Last in table
    expect_equal(isHGNCsymbol(c(NA, "A1CF", NULL, "a2m")), c(FALSE, TRUE, TRUE))
})

test_that("unexpected input does not lead to output", {
    expect_error(isHGNCsymbol(mean), "cannot coerce")
})

# [END]
