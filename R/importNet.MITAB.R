###psuedocode for psi_import###


#A dataframe containing information about whether or not a type of interaction requires the addition of a symmetric edge.
all_interactions = dataframe of all possible interactions and a column called "symmetric", with Boolean values TRUE and FALSE

#Read in all information from the PSI-MITAB file
interactions = read(fname)

#Keep only the unique identifier for A and B, Aliases for A and B, interaction detection methods,
#NCBI tax id for A and B, and confidence scores columns
keep_column = c("Unique Identifier A", "Unique Identifier B", "tax ID A", "tax ID B", "Confidence scores", "interaction_detection")
interactions = interactions[ ,keep_column]

#Remove interactions which dont have both interactors equal to request taxID
interactions[(which(interactions$taxID_A != request_taxID) | which(interactions$taxID_B != request_taxID)), ]

#Remove interactions which have an interaction detection method that is not specified in the list of interaction methods provided by user
interactions[interaction_detection_method %in% method, ]

interactions$unique_identifier_A = FastMap(interactions$unique_identifier_A, "uniprot")
interactions$unique_identifier_B = FastMap(interactions$unique_identifier_B, "uniprot")

#Create normalized confidence scores for InAct interactions
if(database == InAct){
    confidence_scores_InAct = interactions$confidence_scores
    confidence_scores_InAct = gsub("-", 0, confidence_scores)
    normalized_scores_InAct = confidence_scores/max(confidence_scores)
}

#create normalized confidence scores for iRefWeb interactions
if(database == iRefWeb){
    confidence_scores_iRefWeb = interactions$confidence_scores
    for(interaction in confidence_scores_iRefWeb){
        if(interaction$hpr > 50 & interaction$lpr > 50){
            confidence_scores_iRefWeb[interaction] = interaction$np*2
        }
        if(interaction$hpr > 50 & interaction$lpr <= 50){
            confidence_scores_iRefWeb[interaction] = interaction$np*2.5
        }
        if(interaction$hpr <= 50 & interaction$lpr <= 50){
            confidence_scores_iRefWeb[interaction] = interaction$np*3
        }
        if(interaction$hpr == "-"){
            confidence_scores_iRefWeb[interaction] = 0
        }
        if(interaction$lpr == "-"){
            confidence_scores_iRefWeb[interaction] = 0
        }
    }
    normalized_scores_iRefWeb = confidence_scores_iRefWeb/max(confidence_scores_iRefWeb)
}

#create normalized confidence scores for inBioMap interactions
if(database == inBioMap){
    confidence_scores_inBioMap = interactions$confidence_scores(value before |)
    gsub("-", 0, confidence_scores_inBioMap)
    normalized_scores_inBioMap = confidence_scores_inBioMap/max(confidence_scores_inBioMap)
}

#created normalized combined confidence scores
for(unique interaction in InAct, iRefWeb, inBioMap){
    combined_confidence_score[interaction] = (normalized_scores_InAct[interaction] + normalized_scores_iRefWeb[interaction] + normalized_scores_inBioMap[interaction])/3
}

interactions$confidence_scores = combined_confidence_score

#Add symmetric edges
for (interaction in interactions){
    if (symmetric value of interaction in all_interactions == TRUE){
        add_edge(b, a)
    }
}

#Remove duplicated edges and loos
interactions = simplify(interactions)

#remove interactions that dont meet the required cutoff
if(cutofftype == "xS"){
    interactions = interactions[which(interactions$confidence_scores >= val), ]
}
if(cutofftype == "xQ"){
    quantile_value = quantile(interactions$confidence_scores, val)
    interactions = interactions[which(interactions$confidence_scores > quantile_value), ]
}
if(cutofftype == "xN"){
    top = tail(sort(interactions$confidence_scores), val)
    interactions = interactions[which(interactions$scores == top), ]
}

## remove columns which not needed in the graph
interactions = interactions[ , -c("taxID_A", "tax_ID_B", "interaction_detection")]

##rename the columns
colnames(interactions) = c("A", "B", "ConfidenceScore")

##create gG object from interactions
gene_graph = .df2gG(interactions)
