#=============Testing Vectorized Hash Query vs Loop Based Hash Query============

#==================================Design=======================================
# Identified structures/packages for fast ID translation (by pairing of two
# values) are:
#
## 1) Hashed R environments (base, no package required): either by using a loop,
# or by the means of the mget() function, which also has a 'ifnotfound' argument
# that replaces the result vector position with the specified variable (NA in
# this case)
#
## 2) Package 'hash', sets up a hash table as an 'enclosure', and has
# the function 'values()' where it takes a hash package specific hashtable
# and returns the values corresponding to each pair, and generates and error if
# a key does not exist, It also has a 'has.keys()' function that checks for the
# existence of keys in a hashtable, and returns a boolean vector.
#
## 3) Package 'hashmap', creates hashtable and uses standard '[[' bracket
# accessor, automatically replacing elements that are not found with 'NA'.
#
## 4) R lists, it has been observed that are lists or named vectors are faster
# than hash tables for structures with less than 100 elements, which will not be
# the case here as there are far more IDs.

## For benchmarking, the package 'microbenchmark' will be used, which measures
# time on a nanosecond scale (on Macs). It can be given a list of jobs to
# to exectue (using 'list' argument) as well as the number of iterations per
# test sample (using the 'times' argument).
#
## Proposed methods and parameters:
# The samples will be generated using a random string generator function,
# with keys of length 6 and values of length 15. (can be adjusted)
# The size of the tables can be changed, but a size of 25000 is proposed for the
# initial tests.
# All the tables use the same key-value pairs, for consistency.
# Missing keys can also be generated using the same function, and are simply not
# added to the table, but kept aside to be sampled from.
# The sample (query) size can be changed, the query can either be in chunks, or
# all keys at once.
# The percentage of missing keys at each query can also be determined using a
# single variable, across all hashtable tests. The impact of percentage
# missing keys on performance and memory use will be measured.
# The queries (or jobs) are to be generated as unevaluated expressions
# before benchmarking by the use of 'bquote()' function, in order to minimize
# computational costs, and avoid dynamic sampling of keys during the benchmarks.
# Memory use can be benchmarked by the use of R's 'mem_used()', 'mem_changed()'
# as well as the package 'lineprof', which allows the measurement of
# incremental changes in memory.


#=================================Setup=========================================
#
# Install and Load packages
# install.packages("hashmap")
# install.packages("hash")
# install.packages("microbenchmark")
# Idenitified
library(hash)
library(hashmap)
library(microbenchmark)

# Set test parameters
testIterations <- 1000
testSize <- 6000
percentM <- 0.05
keySize <- 6
valueSize <- 15
characters <- c(LETTERS, 0:9)
tableSize <- 25000
set.seed(5643728)

## Create 6 letter keys and 15 letter values to simulate real world scenario
# Random name generator function
keyGenerator <- function(nameSize, characters, count) {
    randomNames <- vector(mode = "character", length = count)
    for (i in 1:count) {
        randomNames[i] <-
            paste(
                sample(characters, size = nameSize, replace = TRUE),
                sep = "",
                collapse = ""
            )
    }
    return(randomNames)
}
# The specific characters do not really matter
valueGenerator <- function(nameSize, characters, count) {
    randomNames <- vector(mode = "character", length = count)
    for (i in 1:count) {
        randomNames[i] <- paste(c("ENSP", sample(
            0:9, size = (nameSize - 4),
            replace = TRUE
        )),
        sep = "",
        collapse = "")
    }
    return(randomNames)
}

## Create random key-value pairs, and misses
Keys <- keyGenerator(keySize, characters, tableSize)
Values <- valueGenerator(valueSize, characters, tableSize)
missKeys <- keyGenerator(keySize, characters, tableSize)


## Set up hash tables and environment

hashPackageTable <- hash::hash(keys = Keys, values = Values)

hashMapPackageTable <-
    hashmap::hashmap(keys = Keys, values = Values)

baseHashTable <- new.env(hash = TRUE, size = tableSize)
for (i in 1:tableSize) {
    baseHashTable[[Keys[i]]] <- Values[i]
}
# env.profile(baseHashTable)

## Different hash tables have different query methods

sampleGenerator <- function() {
    noMiss <- sample(Keys, size = testSize * (1 - percentM))
    miss <- sample(missKeys, size = testSize * percentM)
    return(sample(c(noMiss, miss)))
}

globalQueries <- vector("list", testIterations)
for (i in 1:testIterations) {
    globalQueries[[i]] <- sampleGenerator()
}

# base jobs with misses, note the ifnotfound argument
baseJobs <- vector("list", testIterations)
for (i in 1:testIterations) {
    baseJobs[[i]] <-
        local({
            bquote(mget(.(globalQueries[[i]]), baseHashTable, ifnotfound = NA))
        })
}

# base jobs with loops...
loopFunc <- function(query) {
    result <- vector(mode = "character", length = testSize)
    for (i in 1:testSize) {
        lookup <- baseHashTable[[query[i]]]
        if (is.null(lookup)) {
            out <- NA
        }
        else {
            out <- lookup
        }
        result[i] <- out
    }
}
baseLoopJobs <- vector("list", testIterations)
for (i in 1:testIterations) {
    baseLoopJobs[[i]] <-
        local({
            bquote(loopFunc(.(globalQueries[[i]])))
        })
}

# hashPackage jobs, note that use of for loop was avoided (using lookup and
# checking for null), and instead query keys were checked for existence in the
# hashtable before the actual queries
hashPackageJobs <- vector("list", testIterations)
for (i in 1:testIterations) {
    hashPackageJobs[[i]] <-
        local({
            bquote(globalQueries[[.(i)]][which(has.key(key = globalQueries[[.(i)]], hash = hashPackageTable))])
        })
}



# hashmapPackage jobs, note that it automatically appends NA (does not generate
# an error)
hashmapPackageJobs <- vector("list", testIterations)
for (i in 1:testIterations) {
    hashmapPackageJobs[[i]] <-
        local({
            bquote(hashMapPackageTable[[globalQueries[[.(i)]]]])
        })
}

#==============================Benchmarks=======================================
baseBenchmark <- microbenchmark(list = baseJobs, times = 1)

baseLoopBenchmark <- microbenchmark(list = baseLoopJobs, times = 1)
hashPackageBenchmark <- microbenchmark(list = hashPackageJobs, times = 1)

hashmapPackageBenchmark <- microbenchmark(list = hashmapPackageJobs, times = 1)
rm (list = c(
    "baseJobs",
    "baseLoopJobs",
    "globalQueries",
    "hashmapPackageJobs"
))
#=============================Plots and Results=================================

plot(
    baseBenchmark$time
    ,
    xlab = "Iterations"
    ,
    ylab = "Time (ns)"
    ,
    col = "green"
    ,
    ylim = c(0, 6 * 10 ^ 7)
)
points(hashPackageBenchmark$time
       , col = "blue")
points(hashmapPackageBenchmark$time
       , col = "red")
points(baseLoopBenchmark$time
       , col = "orange")
title(
    main = paste(
        "Comparison of package:hash, package:hashmap & base:env\n",
        "Using values(), [[, mget() and env loop\n",
        "Misses:",
        percentM,
        ", Query size:",
        testSize
    ),
    font.main = 14,
    font = 2
)
legend(
    x = 1,
    y = 6 * 10 ^ 7,
    legend = c("env", "hash", "hashmap", "env loop")
    ,
    cex = 0.4,
    col = c("green", "blue", "red", "orange"),
    pch = 21
)

#===============================Notes===========================================


# [END]
