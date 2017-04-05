### pseudocode for symmetrizeMITAB function

## download PSI controlled vocabulary list of methods from the EMBL-EBI Ontology Lookup Service
## interaction_vocab is the name of the files containing the vocabulary list of interaction methods
all_interaction_methods <- getOntology(interaction_vocab)

## complile list with all experimental interaction detection methods
experimental <- getTermAndAllChildren(all_interaction_methods, "MI:0045")

## complile list with all inference interaction detection methods
inference <- getTermAndAllChildren(all_interaction_methods, "MI:0362")

## complile list with all interaction prediction interaction detection methods
interaction_prediction <- getTermAndAllChildren(all_interaction_methods, "MI:0083")

## is_Symmetric_experimental, is_symmetric_inference_inference, and is_symmetric_interaction_prediction
## are vectors containing boolean values TRUE if the interaction method is symmetric or
## FALSE if the interaction method is not symmetic or NA if there is not enough
## information to determine if the interaction method is symmetric or not.
## These vectors will be compiled by through literacture curation during this task
experimental <- cbind(experimental, is_symmetric_method_experimental)
inference <- cbind(is_symmetric_method_inference)
interaction_prediction <- cbind(is_symmertic_method_interaction_prediction)

## all_interactions is a dataframe cantaining all the interactions, that require checking to see
## if they imply a symmetric interaction. This is the input that would be passed into this function
## from the IMPORT-N modules. user_methods is a vector containing the methods provided by the user.
## Remove interactions which have a method that is not specified in user_methods.
all_interactions[all_interactions$methods %in% user_method, ]

## Symmetric_interactions is a dataframe containing all the interactions in all_interactions
## but with the interactor a and interactor b swapped
symmetric_interactions  = all_interactions
interactor_A = all_interactions$interactor_A
interactor_B = all_interactions$interactor_B
symmetric_interactions$interactor_A = interactor_B
symmetric_interactions$interactor_B = interactor_A

for (i in length(nrows(symmetric_interactions))){
    parents = getAllTermParents(all_interaction_methods, symmetric_interactions$method[i]))
    ## if the interaction method is a type of experimental method, then we search our list of all
    ## experimental methods to find the method
    if (any parents %in% experimental == TRUE){
        x = which(experimental == symmetric_interactions$method[i])
        ## once the method is found we check whether it is symmetric, and if it is not, any interactions
        ## that have this method are removed from the symmetic_interactions dataframe
        if (experimental$is_symmetric_method_experimental[x] != TRUE){
            nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
            symmetric_interactions = symmetric_interactions[nonsymmetric_interactions, ]
        }
    }
    ## if the interaction method is a type of inference method, then we search our list of all
    ## inference methods to find the method
    if (any parents %in% inference == TRUE){
        x = which(inference == symmetric_interactions$method[i])
        ## once the method is found we check whether it is symmetric, and if it is not, any interactions
        ## that have this method are removed from the symmetic_interactions dataframe
        if (inference$is_symmetric_method_experimental[x] != TRUE){
            nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$methos[i])
            symmetric_interactions = symmetric_interactions[nonsymmetric_interactions, ]
        }
    }
    ## if the interaction method is a type of interaction_prediction method, then we search our list of all
    ## interaction_prediction methods to find the method
    if (any parents %in% interaction_prediction == TRUE){
        x = which(interaction_prediction == symmetric_interactions$method[i])
        ## once the method is found we check whether it is symmetric, and if it is not, any interactions
        ## that have this method are removed from the symmetic_interactions dataframe
        if (interaction_prediction$is_symmetric_method_experimental[x] != TRUE){
            nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
            symmetric_interactions = symmetric_interactions[nonsymmetric_interactions, ]
        }
    }
    ## if the method is not experimental, inference, or interaction_prediction method, then
    ## it is an unspecified interaction and we cannot determine if it is symmetric, so we will
    ## remove all interactions with this method
    else{
        nonsymmetric_interactions = which(symmetric_interactions$method == symmetric_interactions$method[i])
        symmetric_interactions = symmetric_interactions[nonsymmetric_interactions, ]
    }
}

##combine and return the list of all_interactions including the added symmetric interactions
all_interactions = rbind(all_interactions, symmetric_interactions)
return all_interactions
