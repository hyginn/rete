symmetrizeMITAB <- function(interactions, methods){

    ## install ontoCAT package if required
    if(require(ontoCAT)){
        print("Trying to install ontoCAT")
        source("https://bioconductor.org/biocLite.R")
        biocLite("ontoCAT")
        library(ontoCAT)
        if(require(ontoCAT)){
            print("ontoCAT could not be installed")
        }
    }
    ##check to see if parameters are valid values
    if(is.null(interactions)){
        errorMessage <- sprintf("No interactions supplied. Please provide interaction data")
        stop(errorMessage)
    }
    if(length(which(is.na(methods))) > 0){
        errorMessage <- sprintf("The interactions supplied contain NAs. Please provide interaction data without NAs")
        stop(errorMessage)
    }
    if(is.data.frame(interactions) == FALSE){
        errorMessage <- sprintf("The interactions parameter supplied is not a dataframe. Please provide interactions dataframe")
        stop(errorMessage)
    }
    if(is.null(methods)){
        errorMessage <- sprintf("No methods supplied. Please provide methods data")
        stop(errorMessage)
    }
    if(length(which(is.na(methods))) > 0){
        errorMessage <- sprintf("Methods contain NAs. Please provide method data without NAs")
        stop(errorMessage)
    }
    if(is.data.frame(methods) == FALSE){
        errorMessage <- sprintf("The methods parameter supplied is not a vector. Please provide methods vector")
        stop(errorMessage)
    }
    ## create interaction method ontology tree
    method_ontology <- getOntology("http://purl.obolibrary.org/obo/mi.owl")
    ## get all the interaction methods which are considered as experimental
    experimental <- getTermAndAllChildrenById(method_ontology, "MI:0045")
    ## get all the interactions which are considered as inference
    inference <- getTermAndAllChildrenById(method_ontology, "MI:0362")
    ## get all the interactions which are considered as interaction prediction
    interaction_prediction <- getTermAndAllChildrenById(method_ontology, "MI:0083")

    ##TO DO: finish manually compiling the is_Symmetric_experimental, is_symmetric_inference_inference, and is_symmetric_interaction_prediction
    is_symmetric_experimental
    is_symmetric_inference
    is_symmetric_interaction_prediction

    ## Add the columns describing if a symmetric interaction is required for the experimental, inference, and interaction_prediction methods
    experimental <- cbind(experimental, is_symmetric_method_experimental)
    inference <- cbind(is_symmetric_method_inference)
    interaction_prediction <- cbind(is_symmertic_method_interaction_prediction)

    ## subset the all the input interactions to include only the interactions with methods specified by the user in input parameter "methods"
    interactions = interactions[interactions$method %in% methods, ]

    ## Create symmetric interactions for all the interactions
    symmetric_interactions  = interactions
    interactor_A = interactions$interactor_A
    interactor_B = interactions$interactor_B
    symmetric_interactions$interactor_A = interactor_B
    symmetric_interactions$interactor_B = interactor_A
    symmetric_interactions$method = interaction$method

    for (i in length(nrows(symmetric_interactions))){
        parents = getAllTermParentsById(method_ontology, symmetric_interactions$method[i])
        ## if the interaction method is a type of experimental method, then we search our list of all
        ## experimental methods to find the method
        if (length(which(parents %in% experimental == TRUE)) > 0){
            x = which(experimental == symmetric_interactions$method[i])
            ## once the method is found we check whether it is symmetric, and if it is not, any interactions
            ## that have this method are removed from the symmetic_interactions dataframe
            if (experimental$is_symmetric_method_experimental[x] != TRUE){
                nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
                symmetric_interactions = symmetric_interactions[-nonsymmetric_interactions, ]
            }
        }
        ## if the interaction method is a type of inference method, then we search our list of all
        ## inference methods to find the method
        if (length(which(parents %in% inference == TRUE)) > 0){
            x = which(inference == symmetric_interactions$method[i])
            ## once the method is found we check whether it is symmetric, and if it is not, any interactions
            ## that have this method are removed from the symmetic_interactions dataframe
            if (inference$is_symmetric_method_experimental[x] != TRUE){
                nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$methos[i])
                symmetric_interactions = symmetric_interactions[-nonsymmetric_interactions, ]
            }
        }
        ## if the interaction method is a type of interaction_prediction method, then we search our list of all
        ## interaction_prediction methods to find the method
        if (length(which(parents %in% interaction_prediction == TRUE)) > 0){
            x = which(interaction_prediction == symmetric_interactions$method[i])
            ## once the method is found we check whether it is symmetric, and if it is not, any interactions
            ## that have this method are removed from the symmetic_interactions dataframe
            if (interaction_prediction$is_symmetric_method_experimental[x] != TRUE){
                nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
                symmetric_interactions = symmetric_interactions[-nonsymmetric_interactions, ]
            }
        }
        ## if the method is not experimental, inference, or interaction_prediction method, then
        ## it is an unspecified interaction and we cannot determine if it is symmetric, so we will
        ## remove all interactions with this method
        else{
            nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
            symmetric_interactions = symmetric_interactions[-nonsymmetric_interactions, ]
        }
    }

interactions = rbind(interactions, symmetric_interactions)
}
