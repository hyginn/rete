# hello.R

# This is an example function which writes 'Hello, world!' to console, and
# returns the first seven Fibonacci numbers.
#

hello <- function() {
  cat("Hello, world!\n")
  return(c(1, 1, 2, 3, 5, 8, 13))
}

# [END]