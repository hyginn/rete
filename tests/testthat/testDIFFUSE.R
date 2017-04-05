#DIFFUSE Tests




# ==== BEGIN SETUP AND PREPARE =================================================
# set up tempdir and tempfiles for output files and to be deleted after tests
OLOG <- as.character(getOption("rete.logfile"))   # save original logfile name
rete::logFileName(fPath = tempdir(), setOption = TRUE)  # make tempdir() the log dir
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}

if (!is.null(getOption("rete.beta"))) {
    storeBeta <- getOption ("rete.beta")
    options("rete.beta" = NULL)
} else {
    storeBeta <- NULL
}
options("rete.beta" = 0.25)
# ==== END SETUP AND PREPARE ===================================================

#Creating Test Data for testDIFFUSE

#========== Make AGG ==============
HGNCsymb<-c("ABC1","DEF2","GHI4","JKL3","LMN5","OPQ6","RST7","UVW8","XYZ9")
VertScores<-c(10,2,4,6,7,8,2,5,9)

aggVertices<-data.frame(HGNC_symbol=HGNCsymb,Gene_Score=VertScores,
                        stringsAsFactors = FALSE)

From<-c("ABC1","ABC1","DEF2","JKL3","DEF2","DEF2","JKL3",
        "GHI4","JKL3","GHI4","GHI4","GHI4","LMN5",
        "XYZ9","LMN5","OPQ6","OPQ6","UVW8","OPQ6",
        "RST7","UVW8","XYZ9")
To<-c("DEF2","JKL3","ABC1","ABC1","JKL3","GHI4","DEF2",
      "DEF2","GHI4","JKL3","LMN5","XYZ9","GHI4",
      "GHI4","OPQ6","LMN5","UVW8","OPQ6",
      "RST7","OPQ6","XYZ9","UVW8")

edgeID<-character(length=length(From))

for (i in 1:length(edgeID)) {
    edgeID[i]<-paste("BNE_",as.character(i),sep="")
}

aggEdges<-data.frame(from=From,to=To,edgeID=edgeID)
#edge data frame created for AGG


AGG<-igraph::graph_from_data_frame(aggEdges, directed = TRUE, vertices = aggVertices)

#setup AGG metadata
attr(AGG, "type") <- "AGG"
attr(AGG, "version") <- "1.0"
attr(AGG, "UUID") <- uuid::UUIDgenerate()

#========= Make matrix W =============
#Creating a reference matrix W to which to compare DIFFUSE output
numVerts<-nrow(aggVertices)

wVect<-numeric(length = numVerts^2)

wMatrixTest<-matrix(nrow = numVerts, ncol = numVerts , data = wVect ,
                    dimnames = list( x = HGNCsymb, y = HGNCsymb ))


#Set up Wij = 1/deg(j) if i interacts with j
wMatrixTest["DEF2","ABC1"] <- 1/2
wMatrixTest["JKL3","ABC1"] <- 1/2
wMatrixTest[c(1,3,4),"DEF2"] <- 1/3
wMatrixTest[c(2,4,5,9),"GHI4"] <- 1/4
wMatrixTest[c(1,2,3),"JKL3"]<-1/3
wMatrixTest[c(3,6),"LMN5"] <-1/2
wMatrixTest[c(5,7,8),"OPQ6"]<-1/3
wMatrixTest[6,"RST7"]<-1
wMatrixTest[c(6,9),"UVW8"]<-1/2
wMatrixTest[c(3,8),"XYZ9"] <- 1/2

#======= Perform subsequent steps to make matrix F and matrix E =========

idMat<-diag(nrow = numVerts , ncol = numVerts)

Beta<-0.25

preInversion <- idMat - (1-Beta)*wMatrixTest

#Perform inversion, multiply by Beta to obtain F (matF)

matFRef<-Beta*solve(preInversion)

#Create a reference matrix E (post multiply matFRef by heat scores)

heatVectRef<-aggVertices$Gene_Score

heatDiagRef<-diag(heatVectRef)

matERef<-matFRef%*%heatDiagRef

aggEdges[,"Influence"] <- numeric(length = nrow(aggEdges))

#======== Assign influence to edges for eggRef (reference EGG) =========

#Manually Assigning Influence

aggEdges$Influence[1] <- matERef[2,1] #ABC1 to DEF2
aggEdges$Influence[2] <- matERef[4,1] #ABC1 to JKL3
aggEdges$Influence[3] <- matERef[1,2] #DEF2 to ABC1
aggEdges$Influence[4] <- matERef[1,4] #JKL3 to ABC1
aggEdges$Influence[5] <- matERef[4,2] #DEF2 to JKL3
aggEdges$Influence[6] <- matERef[3,2] #DEF2 to GHI4
aggEdges$Influence[7] <- matERef[2,4] #JKL3 to DEF2
aggEdges$Influence[8] <- matERef[2,3] #GHI4 to DEF2
aggEdges$Influence[9] <- matERef[3,4] #JKL3 to GHI4
aggEdges$Influence[10] <- matERef[4,3] #GHI4 to JKL3
aggEdges$Influence[11] <- matERef[5,3] #GHI4 to LMN5
aggEdges$Influence[12] <- matERef[9,3] #GHI4 to XYZ9
aggEdges$Influence[13] <- matERef[3,5] #LMN5 to GHI4
aggEdges$Influence[14] <- matERef[3,9] #XYZ9 to GHI4
aggEdges$Influence[15] <- matERef[6,5] #LMN5 to OPQ6
aggEdges$Influence[16] <- matERef[5,6] #OPQ6 to LMN5
aggEdges$Influence[17] <- matERef[8,6] #OPQ6 to UVW8
aggEdges$Influence[18] <- matERef[6,8] #UVW8 to OPQ6
aggEdges$Influence[19] <- matERef[7,6] #OPQ6 to RST7
aggEdges$Influence[20] <- matERef[6,7] #RST7 to OPQ6
aggEdges$Influence[21] <- matERef[9,8] #UVW8 to XYZ9
aggEdges$Influence[22] <- matERef[8,9] #XYZ9 to UVW8

#======= Make eggRef (reference EGG), attach metadata ==========

eggRef <- igraph::graph_from_data_frame(aggEdges, directed = TRUE, vertices= aggVertices)

#======================= Conduct Tests (Finally) ============================

testthat::context ("Check that errors are thrown with erroneous input")

testthat::test_that("Erroneous input -> Error message", {

    #=======Make bad AGGs to see if erroneous inputs induce errors =========

    #Here, we make a whole bunch of bad AGGs that should induce errors

    badAGG1<-AGG #Remove 1st 4 edges; throw error because unused vertices
    badAGG2<-AGG #Remove 1st 4 vertices; throw error because unused edges.
    #Turns out that if you remove vertices as done for badAGG2, igraph will
    #remove all edges including them
    badAGG3<-AGG #Remove Scores
    badAGG4<-AGG #Remove HGNC symbols from vertices
    badAGG5<-AGG #Insert character vector into one of the 'scores'
    badAGG6<-AGG #Insert logical TRUE into one of the 'scores'
    badAGG7<-AGG #Insert integer "1" for one of the 'scores'
    badAGG8<-AGG #Insert factor "1.5" for one of the 'scores'
    badAGG9<-AGG #Give an EGG version number and other metadata. Don't want somebody mistakenly running EGG
    #through pipelin
    badAGG10<-AGG #Give AGG Edge score attribute. Don't want to run EGG through DIFFUSE again

    badAGG1 <- igraph::delete_edges(badAGG1, 1:4)
    badAGG2 <- igraph::delete_vertices(badAGG2, 1:4)
    badAGG3 <- igraph::remove.vertex.attribute(badAGG3,"Gene_Score")
    badAGG4 <- igraph::remove.vertex.attribute(badAGG4,"name")
    badAGG5 <- igraph::vertex_attr(badAGG5,"Gene_Score", igraph::V(badAGG5))[3] <-"cow"
    badAGG6 <- igraph::vertex_attr(badAGG6,"Gene_Score", igraph::V(badAGG6))[3] <-TRUE
    badAGG7 <- igraph::vertex_attr(badAGG7,"Gene_Score", igraph::V(badAGG7))[3] <-as.integer(1)
    badAGG8 <- igraph::vertex_attr(badAGG8,"Gene_Score", igraph::V(badAGG8))[3] <-as.factor(1.5)

    attr(badAGG9, "type") <- "EGG"

    badAGG10 <- igraph::vertex_attr(badAGG10,"Influence")<-numeric(length=length(igraph::V(badAGG10)))






    #expect an error thrown if:
    testthat::expect_error(DIFFUSE(AGG, algorithm="Lies" , param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    #algorithm character vector not recognized
    testthat::expect_error(DIFFUSE(AGG, algorithm=1 , param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    #algorithm is numeric
    testthat::expect_error(DIFFUSE(AGG, algorithm=TRUE , param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    #algorithm is a logical
    testthat::expect_error(DIFFUSE(AGG, algorithm=as.factor("Leis") ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    #algorithm is a factor
    testthat::expect_error(DIFFUSE(AGG, algorithm = "Leis" ,
                         param <- list(getOption("rete.alpha")),
                         silent = TRUE, writeLog = TRUE))
    #option does not exist for param
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list("puppies"),
                         silent = TRUE, writeLog = TRUE))
    #value in param not recognized
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" , param <- 0.25,
                         silent = TRUE, writeLog = TRUE))
    #param not in list format
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = "FALSE", writeLog = TRUE))
    #silent is a character vector

    #testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
     #                    param <- list(getOption("rete.beta")),
      #                   silent = as.factor(FALSE), writeLog = FALSE))
    #silent is a factor: test removed
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = 2, writeLog = TRUE))
    #silent is a double numeric
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = "TRUE"))
    #writeLog is a character vector
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = as.factor(TRUE)))
    #writeLog is a factor
    testthat::expect_error(DIFFUSE(AGG, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = 2))
    #writeLog is a double numeric
    testthat::expect_error(DIFFUSE(badAGG1, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    #see lines11-21 for description of bad AGG inputs that should cause
    #an error
    #testthat::expect_error(DIFFUSE(badAGG2, algorithm="Leis" ,
                         #param <- list(getOption("rete.beta")),
                         #silent = FALSE, writeLog = TRUE))
    #it turns out that iGraph removes edges if you remove their vertices,
    #rendering a test to see what happens when vertices removed unnecessary

    testthat::expect_error(DIFFUSE(badAGG3, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG4, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG5, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG6, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG7, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG8, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG9, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))
    testthat::expect_error(DIFFUSE(badAGG10, algorithm="Leis" ,
                         param <- list(getOption("rete.beta")),
                         silent = TRUE, writeLog = TRUE))

})

testthat::context("Check that edge scores assigned to EGG produced in test match
        those in manually constructed reference")

testthat::test_that("EGG made correctly", {
    EGG<-DIFFUSE(AGG, algorithm = "Leis",
                 param = list(getOption("rete.beta")),
                 silent = TRUE, writeLog = TRUE)
    eggEdges <- igraph::as_data_frame(EGG, what = "edges")
    eggRefEdges<-igraph::as_data_frame(eggRef, what = "edges")

    fromLogiVect <- eggEdges$from == eggRefEdges$from
    toLogiVect <- eggEdges$to == eggRefEdges$to
    edgeIDLogiVect <- eggEdges$edgeID == eggRefEdges$edgeID
    inflLogiVect <- eggEdges$Influence == eggRefEdges$Influence

    testthat::expect_equal(sum(as.numeric(as.logical(fromLogiVect))), length(fromLogiVect))
    testthat::expect_equal(sum(as.numeric(as.logical(toLogiVect))), length(toLogiVect))
    testthat::expect_equal(sum(as.numeric(as.logical(edgeIDLogiVect))), length(edgeIDLogiVect))
    testthat::expect_equal(sum(as.numeric(as.logical(inflLogiVect))), length(inflLogiVect))
    # for each comparison of columns, we expect all rows to have value of true


    #Have we got the log file?

    # correct log file
    thisLog <- readLines(logName)
    #testthat::expect_true(grepl("event \\| title  \\| DIFFUSE",     thisLog[1]))


    # call
    #testthat::expect_true(grepl("^event \\| call   \\| DIFFUSE\\(", thisLog[3]))
    #testthat::expect_true(grepl("algorithm = \"Leis\"",                    thisLog[3]))
    #testthat::expect_true(grepl("AGG = \"AGG\"",                       thisLog[3]))
    #testthat::expect_true(grepl("param = [[1]][1] 0.25", thisLog[3]))
    #testthat::expect_true(grepl("silent = TRUE",                             thisLog[3]))
    #testthat::expect_true(grepl("writeLog = TRUE",                           thisLog[3]))
    #testthat::expect_true(grepl(")$",                                        thisLog[3]))

    #notes


})

#To do: have a way of checking that metadata is assigned correctly

testthat::test_that("Metadata assigned correctly", {
    EGG<-DIFFUSE(AGG, algorithm = "Leis",
                 param <- list(getOption("rete.beta")),
                 silent = TRUE, writeLog = TRUE)
    eggEdges <- igraph::as_data_frame(EGG, what = "edges")
    eggRefEdges<-igraph::as_data_frame(eggRef, what = "edges")

    testthat::expect_equal(attr(EGG,"type"), "EGG")
    testthat::expect_equal(attr(EGG,"version"), "1.0")
    testthat::expect_equal(length(.checkArgs(attr(EGG, "UUID"),
                                      uuid::UUIDgenerate() ) ),
                           0)


})

#Below test was totally plagiarized from diffuseV3.R, but better not to
#re-engineer the wheel
testthat::test_that("DIFFUSE doesn't leak output if silent = TRUE", {
    testF <- tempfile()
    capture.output(dummy <- DIFFUSE(AGG = AGG, algorithm = "Leis",
                   param = list(getOption("rete.beta")),
                   silent = TRUE,
                   writeLog = FALSE),
                   file = testF)
    testthat::expect_equal(length(readLines(testF)), 0)
    testthat::expect_true(file.remove(testF))
})

#To do:
#Come up with a test to see that nothing is written to the log if writeLog is TRUE

# ==== BEGIN TEARDOWN AND RESTORE ==============================================
logName <- unlist(getOption("rete.logfile"))
if (file.exists(logName)) { file.remove(logName)}

if ( !(is.null(storeBeta))) {
    options("rete.beta" = storeBeta)
} else {
    options("rete.beta"= NULL)
}

options("rete.logfile" = OLOG)
# ==== END  TEARDOWN AND RESTORE ===============================================

#[END]
