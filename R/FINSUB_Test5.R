#Test 5

#================Test5===============================
#Here, we're just checking that all the checking conditions work. See the 1st 50 or so lines of code in
#findsub.R

folder<-"C:/Users/HPDM1/Documents/CanadaUofT/4thYear/BCB420/ekplektoR/R/"
fPath<-paste(folder,"FINDSUB_Test_Data.R",sep="")
source(fPath) #executing script generating test data. See script for variables generated, network structure

context("Checking the checks")

test_that("Checking the checks", {
                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                       noLog = TRUE,silent=NULL)) #silent is null
                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                        noLog = TRUE,silent=2)) #silent is numeric
                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=1,
                                       noLog = TRUE,silent="TRUE")) #silent with character

                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                        noLog = NULL,silent=FALSE)) #noLog is null
                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                   noLog = "NULL",silent=FALSE)) #noLog is Character
                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                   noLog = 43,silent=FALSE)) #noLog is numeric

                expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=-3,
                                                   noLog = TRUE,silent=TRUE))
                expect_error(outputGraphs<-findsub(method="Blargh",EGG,minOrd=3, #character method
                                                        noLog = TRUE,silent=TRUE))

                expect_error(outputGraphs<-findsub(method=NULL,EGG,minOrd=3, #null method
                                                                    noLog = TRUE,silent=TRUE))

                 expect_error(outputGraphs<-findsub(method=2,EGG,minOrd=3, #numeric method
                                                        noLog = TRUE,silent=TRUE))

               expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=NULL, #null minord
                                                      noLog = TRUE,silent=TRUE))



                      expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=2.2,
                                                             noLog = TRUE,silent=TRUE)) #minOrd with non-integer

       expect_error(outputGraphs<-findsub(method=2,EGG,minOrd="Pizza",
                                              noLog = TRUE,silent=TRUE)) #MinOrd w character
       expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=1,
                                          noLog = TRUE,silent=TRUE)) #minOrd<2

})

#An added test: checking the version number

context("Checking that you can check version number")

igraph::graph_attr(EGG,"EGGversion")<-"1.2"

testthat::test_that("Wrong version number -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                                    noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"EGGversion")<-1

testthat::test_that("EGGversion of type double -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"EGGversion")<-as.integer(1)

testthat::test_that("EGGversion of type integer -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::remove.graph.attribute(EGG,"EGGversion")

testthat::test_that("EGGversion nonexistent -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

context("Have a valid delta value")

igraph::graph_attr(EGG,"delta")<-"puppies"

testthat::test_that("Delta a character -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::graph_attr(EGG,"delta")<-as.integer(1)

testthat::test_that("Delta an integer -> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})

igraph::remove.graph.attribute(EGG,"delta")

testthat::test_that("No Delta-> Error message", {
    expect_error(outputGraphs<-findsub(method="Leis",EGG,minOrd=3,
                                       noLog = TRUE,silent=TRUE))
})