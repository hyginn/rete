#Test 5

#================Test5===============================
#Here, we're just checking that all the checking conditions work. See the 1st 50 or so lines of code in
#FINDSUBv1_2.R

folder<-"C:/Users/HPDM1/Documents/CanadaUofT/4thYear/BCB420/ekplektoR/R/"
fPath<-paste(folder,"FINDSUB_Test_Data.R",sep="")
source(fPath) #executing script generating test data. See script for variables generated, network structure

context("Checking the checks")

test_that("Checking the checks", {
                expect_error(outputGraphs<-FINDSUBv1_2(MethType="Leis",Thresh=4,inputGraph,minOrd=3,
                                                       logResults = TRUE,silent=NULL))

                expect_error(outputGraphs<-FINDSUBv1_2(MethType="Leis",Thresh=4,inputGraph,minOrd=3,
                                                        logResults = NULL,silent=FALSE))

                expect_error(outputGraphs<-FINDSUBv1_2(MethType="Leis",Thresh=4,inputGraph,minOrd=1,
                                                       logResults = TRUE,silent=TRUE))
                expect_error(outputGraphs<-FINDSUBv1_2(MethType="Blargh",Thresh=4,inputGraph,minOrd=3,
                                                        logResults = TRUE,silent=TRUE))
                expect_error(outputGraphs<-FINDSUBv1_2(MethType=NULL,Thresh=4,inputGraph,minOrd=3,
                                                                    logResults = TRUE,silent=TRUE))
                 expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=4,inputGraph,minOrd=3,
                                                        logResults = TRUE,silent=TRUE))
              expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=4,inputGraph,minOrd=3,
                                                                     logResults = TRUE,silent=TRUE))
               expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=-1,inputGraph,minOrd=3,
                                                      logResults = TRUE,silent=TRUE))
            expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh="Barn",inputGraph,minOrd=3,
                                                   logResults = TRUE,silent=TRUE))

         expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=2,inputGraph,minOrd=1,
                                                logResults = TRUE,silent=TRUE))

                      expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=2,inputGraph,minOrd=2.2,
                                                             logResults = TRUE,silent=TRUE))
       expect_error(outputGraphs<-FINDSUBv1_2(MethType=2,Thresh=2,inputGraph,minOrd="Pizza",
                                              logResults = TRUE,silent=TRUE))
})


